include { FASTQ_ALIGN_BWA                               } from '../../nf-core/fastq_align_bwa/main'
include { GATK4_REMOVEDUPLICATES                        } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_SORT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGNMENTS   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                             } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS                                } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUPLICATED } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_DEDUP_BWAMEM {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_bwamem_index      // channel: [ val(meta), [ bwam index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments
    use_gpu              // boolean: whether to use GPU or CPU for bwamem alignment

    main:

    ch_alignment                     = Channel.empty()
    ch_alignment_index               = Channel.empty()
    ch_flagstat                      = Channel.empty()
    ch_stats                         = Channel.empty()
    ch_picard_metrics                = Channel.empty()
    ch_multiqc_files                 = Channel.empty()
    ch_versions                      = Channel.empty()

    FASTQ_ALIGN_BWA (
        ch_reads,
        ch_bwamem_index,
        ch_fasta,
        true, // val_sort_bam hardcoded to true
        false // use_gpu hardcoded to false
    )
    ch_alignment        = FASTQ_ALIGN_BWA.out.bam
    ch_alignment_index  = FASTQ_ALIGN_BWA.out.bai
    ch_stats            = FASTQ_ALIGN_BWA.out.stats    // channel: [ val(meta), path(stats) ]
    ch_flagstat         = FASTQ_ALIGN_BWA.out.flagstat // channel: [ val(meta), path(flagstat) ]
    ch_idxstats         = FASTQ_ALIGN_BWA.out.idxstats // channel: [ val(meta), path(idxstats) ]
    ch_versions         = FASTQ_ALIGN_BWA.out.versions

    if (!skip_deduplication) {
        /*
         * Run GATK4 (picard) MarkDuplicates with the --REMOVE_DUPLICATES true flag
         */
        GATK4_REMOVEDUPLICATES (
            ch_alignment
        )

        /*
         * Run samtools index on deduplicated alignment
         */
        SAMTOOLS_INDEX_DEDUPLICATED (
            GATK4_REMOVEDUPLICATES.out.bam
        )
        ch_alignment       = GATK4_REMOVEDUPLICATES.out.bam
        ch_alignment_index = SAMTOOLS_INDEX_DEDUPLICATED.out.bai
        ch_picard_metrics  = GATK4_REMOVEDUPLICATES.out.metrics
        ch_versions        = ch_versions.mix(GATK4_REMOVEDUPLICATES.out.versions)
        ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_DEDUPLICATED.out.versions)

        // Do we want to run samtools stats on the deduplicated alignment?
    }

    /*
     * Collect MultiQC inputs
     */
    ch_multiqc_files = ch_picard_metrics.collect{ meta, metrics -> metrics }
                        .mix(ch_flagstat.collect{ meta, flagstat -> flagstat })
                        .mix(ch_stats.collect{ meta, stats -> stats  })
                        .mix(ch_idxstats.collect{ meta, stats -> stats  })

    emit:
    bam                           = ch_alignment                     // channel: [ val(meta), [ bam ]       ]
    bai                           = ch_alignment_index               // channel: [ val(meta), [ bai ]       ]
    samtools_flagstat             = ch_flagstat                      // channel: [ val(meta), [ flagstat ]  ]
    samtools_stats                = ch_stats                         // channel: [ val(meta), [ stats ]     ]
    samtools_index_stats          = ch_idxstats                      // channel: [ val(meta), [ idxstats ]  ]
    picard_metrics                = ch_picard_metrics                // channel: [ val(meta), [ metrics ]   ]
    multiqc                       = ch_multiqc_files                 // channel: [ *{html,txt}              ]
    versions                      = ch_versions                      // channel: [ versions.yml             ]
}
