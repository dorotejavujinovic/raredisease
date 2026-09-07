//
// A variant caller workflow for parabricks deepvariant
//

include { ADD_VARCALLER_TO_BED                    } from '../../../modules/local/add_varcallername_to_bed'
include { BCFTOOLS_ANNOTATE                       } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_PB   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_PB } from '../../../modules/nf-core/bcftools/norm/main'
include { GLNEXUS                                 } from '../../../modules/nf-core/glnexus/main'
include { PARABRICKS_DEEPVARIANT                  } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { TABIX_BGZIP                             } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_PARABRICKS         } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_SNV_DEEPVARIANT_PARABRICKS {
    take:
        ch_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_case_info       // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header  // channel: [mandatory] [ path(header) ]
        ch_genome_chrsizes // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_par_bed                   // channel: [optional] [ val(meta), path(bed) ] - unused: pbrun handles haploid contigs via ext.args
        ch_target_bed                // channel: [mandatory] [ val(meta), path(bed), path(index) ]
        val_analysis_type            // string:  'wgs' or 'wes'
        val_skip_split_multiallelics // boolean

    main:

        ch_publish = channel.empty()

        if (val_analysis_type.equals("wes")) {
            TABIX_BGZIP(ch_target_bed.map{meta, gzbed, _index -> return [meta, gzbed]})
            ch_bam_bai
                .combine (TABIX_BGZIP.out.output.map {_meta, bed -> return bed})
                .set { ch_parabricks_in }
        } else if (val_analysis_type.equals("wgs")) {
            ch_bam_bai
                .map { meta, bam, bai ->
                        return [meta, bam, bai, []] }
                .set { ch_parabricks_in }
        }

        PARABRICKS_DEEPVARIANT ( ch_parabricks_in, ch_genome_fasta )

        // pbrun writes the gvcf without an index; downstream consumers of
        // gvcf_tabix expect one (native deepvariant emits gvcf_tbi itself).
        TABIX_PARABRICKS ( PARABRICKS_DEEPVARIANT.out.gvcf )

        PARABRICKS_DEEPVARIANT.out.gvcf
            .map{ _meta, gvcf -> gvcf}
            .toSortedList{a, b -> a.name <=> b.name}
            .toList()
            .set { ch_file_list }

        ch_case_info
            .combine(ch_file_list)
            .map {meta, gvcf -> return [meta, gvcf, []]}
            .set { ch_gvcfs }

        GLNEXUS ( ch_gvcfs, [[:],[]] )

        ch_split_multi_in = GLNEXUS.out.bcf
                            .map{ meta, bcf ->
                                    return [meta, bcf, []] }

        if (!val_skip_split_multiallelics) {
            SPLIT_MULTIALLELICS_PB (ch_split_multi_in, ch_genome_fasta)
            ch_remove_dup_in = SPLIT_MULTIALLELICS_PB.out.vcf
                                .map{ meta, vcf ->
                                        return [meta, vcf, []] }
        } else {
            ch_remove_dup_in = ch_split_multi_in
        }
        REMOVE_DUPLICATES_PB (ch_remove_dup_in, ch_genome_fasta)

        ch_genome_chrsizes.flatten().map{chromsizes ->
            return [[id:'parabricks_deepvariant'], chromsizes]
            }
            .set { ch_varcallerinfo }

        ADD_VARCALLER_TO_BED (ch_varcallerinfo).gz_tbi
            .map{_meta,bed,tbi -> return [bed, tbi]}
            .set{ch_varcallerbed}

        REMOVE_DUPLICATES_PB.out.vcf
            .join(REMOVE_DUPLICATES_PB.out.tbi)
            .combine(ch_varcallerbed)
            .combine(ch_foundin_header)
            .map { meta, vcf, vcf_tbi, bed, bed_tbi, hdr -> return [meta, vcf, vcf_tbi, bed, bed_tbi, [], hdr, []] }
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in)

    emit:
        gvcf       = PARABRICKS_DEEPVARIANT.out.gvcf // channel: [ val(meta), path(gvcf)] - already compressed
        gvcf_tabix = TABIX_PARABRICKS.out.index      // channel: [ val(meta), path(gvcf_tbi)] - tabix/tabix emits 'index'
        publish    = ch_publish                      // channel: [ val(meta), path(report) ] - pbrun produces no vcf_stats_report
        tabix      = BCFTOOLS_ANNOTATE.out.tbi       // channel: [ val(meta), path(tbi) ]
        vcf        = BCFTOOLS_ANNOTATE.out.vcf       // channel: [ val(meta), path(vcf) ]
}
