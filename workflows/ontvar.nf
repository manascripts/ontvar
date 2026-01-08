/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ontvar_pipeline'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { SNIFFLES } from '../modules/nf-core/sniffles/main'
include { CUTESV   } from '../modules/nf-core/cutesv/main'
include { SEVERUS as SEVERUS_WITH_CONTROL } from '../modules/nf-core/severus/main'
include { SEVERUS as SEVERUS_NO_CONTROL   } from '../modules/nf-core/severus/main'
include { RENAME_VCF } from '../modules/local/rename_vcf/main'
include { RENAME_VCF_HEADERS as RENAME_VCF_HEADERS_SNIFFLES } from '../modules/local/rename_vcf_headers/main'
include { RENAME_VCF_HEADERS as RENAME_VCF_HEADERS_CUTESV   } from '../modules/local/rename_vcf_headers/main'
include { RENAME_VCF_HEADERS as RENAME_VCF_HEADERS_SEVERUS  } from '../modules/local/rename_vcf_headers/main'
include { JASMINESV as JASMINESV_SAMPLE } from '../modules/nf-core/jasminesv/main'
include { JASMINE_HEADER_FIX } from '../modules/local/jasmine_header_fix/main'
include { JASMINESV as JASMINESV_COHORT } from '../modules/nf-core/jasminesv/main'
include { FILTER_CHR } from '../modules/local/filter_chr/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_COHORT_RAW } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_COHORT    } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_PER_SAMPLE_RAW    } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_PER_SAMPLE    } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../modules/nf-core/annotsv/installannotations/main'
include { UNTAR as UNTAR_ANNOTSV } from '../modules/nf-core/untar/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLERS          } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLER_MERGED    } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLER_MERGED_FILTERED  } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_ANNOTATED } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_FILTERED  } from '../modules/local/summarize_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_RAW_CALLERS          } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_CONSENSUS            } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_FILTERED             } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_COHORT_ANNOTATED     } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_COHORT_FILTERED      } from '../modules/local/plot_sv_counts/main'
include { SVDB_QUERY as SVDB_QUERY_SAMPLE } from '../modules/nf-core/svdb/query/main'
include { SVDB_QUERY as SVDB_QUERY_COHORT } from '../modules/nf-core/svdb/query/main'
include { BCFTOOLS_VIEW as CALLER_SUPPORT_FILTER } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as AF_FILTER } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as AF_FILTER_COHORT } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_SORT as SORT_VCF } from '../modules/nf-core/bcftools/sort/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONTVAR {

    take:
        ch_samplesheet // channel: samplesheet read in from --input
        ch_output_dir // channel: output directory from --outdir
        reference
        annotsv_annotations

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_sample_info = ch_samplesheet
    // ch_sample_info now contains: [group_id, sample_id, sample_type, input_type, input_path]
    //                               [   0   ,    1    ,     2     ,    3      ,     4     ]
    
    // Separate by input type
    fastq_samples = ch_sample_info.filter { it[3] == 'fastq' }
    bam_samples = ch_sample_info.filter { it[3] == 'bam' }

    // ──────────────────────────────────────────────────────────────────────
    // FASTQ CONCATENATION (if multiple files per sample)
    // ──────────────────────────────────────────────────────────────────────
    // Process directory inputs - collect FASTQ files
    fastq_samples_dir = fastq_samples
        .filter { it[4] && file(it[4]).isDirectory() }
        .map { it ->
            def group_id = it[0]
            def sample_id = it[1]
            def fastq_path = it[4]
            def pattern = "${fastq_path}/*.{fastq,fq,fastq.gz,fq.gz}"
            def fastqs = file(pattern)
            tuple([id: "${group_id}", sample: sample_id], fastqs)
        }
    
    // Process single file inputs
    fastq_samples_file = fastq_samples
        .filter { it[4] && !file(it[4]).isDirectory() }
        .map { it ->
            def group_id = it[0]
            def sample_id = it[1]
            def fastq_path = it[4]
            tuple([id: "${group_id}", sample: sample_id], file(fastq_path))
        }

    // Concatenate FASTQs from directories
    CAT_FASTQ(fastq_samples_dir)

    // Merge concatenated and single-file FASTQs for minimap2
    minimap2_input = CAT_FASTQ.out.reads
        .mix(fastq_samples_file)
    
    // ──────────────────────────────────────────────────────────────────────
    // ALIGNMENT (minimap2 for long-read sequencing)
    // ──────────────────────────────────────────────────────────────────────

    MINIMAP2_ALIGN(
        minimap2_input,                                                          // Input 1: [meta, reads]
        Channel.value(tuple([id: "reference"], file(params.reference))),        // Input 2: [meta2, reference]
        Channel.value(true),                                                     // Input 3: bam_format (true for BAM output)
        Channel.value('bai'),                                                    // Input 4: bam_index_extension ('bai' or 'csi')
        Channel.value(false),                                                    // Input 5: cigar_paf_format (false, not needed for BAM)
        Channel.value(true) 
    )

    // Merge aligned and original BAMs with sample_id as key
    aligned_bams = MINIMAP2_ALIGN.out.bam
        .map { meta, bam ->
            def sample_id = meta.id
            tuple(sample_id, bam)
        }
    
    original_bams = bam_samples
        .map { it ->
            def sample_id = it[1]  // sample_id is index 1
            def bam_path = it[4]
            tuple(sample_id, file(bam_path))
        }

    // Merge aligned and original BAMs (keyed by sample_id)
    all_bams = aligned_bams.mix(original_bams)

    // Map cases and controls with sample_id, then join
    cases_with_bams = ch_sample_info
        .filter { it[2] == 'case' }
        .map { it -> tuple(it[1], it[0]) }  // [sample_id, group_id]
        .join(all_bams, by: 0)               // [sample_id, group_id, bam]
        .map { group_id, bam ->
            tuple(group_id, bam)             // [group_id, bam]
        }

    controls_with_bams = ch_sample_info
        .filter { it[2] == 'control' }
        .map { it -> tuple(it[1], it[0]) }  // [sample_id, group_id]
        .join(all_bams, by: 0)               // [sample_id, group_id, bam]
        .map { group_id, bam ->
            tuple(group_id, bam)             // [group_id, bam]
        }

    sv_input = cases_with_bams
        .join(controls_with_bams, by: 0, remainder: true)
        .map { group_id, case_bam, control_bam ->
            tuple(group_id, case_bam, control_bam ?: null)
        }

    // ──────────────────────────────────────────────────────────────────────
    // SV Calling
    // ──────────────────────────────────────────────────────────────────────

    // SNIFFLES
    sniffles_input = sv_input
        .map { it ->
            def group_id = it[0]
            tuple([id: "${group_id}", sample: group_id], it[1], file("${it[1]}.bai"))
        }

    SNIFFLES(
        sniffles_input,                                                         // Input 1: [meta, bam, bai]
        Channel.value(tuple([id: "reference"], file(reference))),               // Input 2: [meta, fasta]
        Channel.value(tuple([id: "tandem"], file(params.tandem_repeats))),      // Input 3: [meta, tandem_file]
        Channel.value(true),                                                    // Input 4: vcf_output
        Channel.value(false)                                                    // Input 5: snf_output
    )

    // CUTESV
    cutesv_input = sv_input
        .map { it ->
            def group_id = it[0]
            tuple([id: "${group_id}", sample: group_id], it[1], file("${it[1]}.bai"))
        }

    CUTESV(
        cutesv_input,
        Channel.value(tuple([id: "reference"], file(reference)))
    )

    // SEVERUS
    severus_with_control_input = sv_input.filter { it[2] }
        .map { it ->
            def group_id = it[0]
            def case_bam = it[1]
            def control_bam = it[2]
            tuple([id: "${group_id}", sample: "${group_id}", has_control: true, control_bam: control_bam], case_bam, file("${case_bam}.bai"), control_bam, file("${control_bam}.bai"), [])
        }

    severus_no_control_input = sv_input.filter { !it[2] }
        .map { it ->
            def group_id = it[0]
            def case_bam = it[1]
            tuple([id: "${group_id}", sample: "${group_id}", has_control: false], case_bam, file("${case_bam}.bai"), [], [], [])
        }

    SEVERUS_WITH_CONTROL(
        severus_with_control_input,                                             // Input 1: [meta, target_bam, target_bai, control_bam, control_bai, vcf]
        Channel.value(tuple([id: "vntr"], file(params.vntr_bed)))               // Input 2: [meta, vntr_bed]
    )

    SEVERUS_NO_CONTROL(
        severus_no_control_input,                                               // Input 1: [meta, target_bam, target_bai, control_bam, control_bai, vcf]
        Channel.value(tuple([id: "vntr"], file(params.vntr_bed)))               // Input 2: [meta, vntr_bed]
    )

    severus_vcfs = SEVERUS_WITH_CONTROL.out.somatic_vcf
        .mix(SEVERUS_NO_CONTROL.out.somatic_vcf)
        .map { meta, vcf ->
            tuple(meta, vcf, 'severus')
        } | RENAME_VCF

    // ──────────────────────────────────────────────────────────────────────
    // FIX SAMPLE NAMES in VCF HEADERS
    // ──────────────────────────────────────────────────────────────────────

    sniffles_renamed_vcfs = SNIFFLES.out.vcf
        .map { meta, vcf -> tuple(meta, vcf) } | RENAME_VCF_HEADERS_SNIFFLES

    cutesv_renamed_vcfs = CUTESV.out.vcf
        .map { meta, vcf -> tuple(meta, vcf) } | RENAME_VCF_HEADERS_CUTESV

    severus_renamed_vcfs = severus_vcfs
        .map { meta, vcf -> tuple(meta, vcf) } | RENAME_VCF_HEADERS_SEVERUS

    all_caller_vcfs = sniffles_renamed_vcfs.mix(cutesv_renamed_vcfs, severus_renamed_vcfs)
        .map { meta, vcf -> tuple(meta, vcf) }

    // Raw caller summaries
    all_caller_vcfs_for_summary = all_caller_vcfs
        .map { meta, vcf -> vcf }
        .collect()
        .map { vcf_list -> tuple([id: "raw_caller_summary"], vcf_list) }

    SUMMARIZE_CALLERS(
        all_caller_vcfs_for_summary,
        Channel.value("raw_calls")
    )

    PLOT_RAW_CALLERS(
        SUMMARIZE_CALLERS.out.json
            .map { meta, json -> tuple([id: "raw_callers_plot"], [json]) },
        Channel.value("Raw Caller SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // Gather SV caller outputs per sample
    // ──────────────────────────────────────────────────────────────────────

    sv_calls_by_sample = all_caller_vcfs
        .map { meta, vcf ->
            def group_id = meta.sample ?: meta.id  // Use sample field first, fallback to id
            tuple(group_id, vcf)
        }
        .groupTuple(by: 0)

    // ──────────────────────────────────────────────────────────────────────
    // Run Jasmine to merge SVs from callers per sample
    // ──────────────────────────────────────────────────────────────────────

    jasminesv_sample_input = sv_calls_by_sample
        .filter { meta, vcf_list -> vcf_list.size() > 0 }
        .map { meta, vcf_list ->
            def group_id = meta
            tuple([id: group_id, sample: group_id, step: "consensus"], vcf_list, [], [])
        }

    // Prepare Jasmine input channels (per-sample)
    ch_jasmine_sample_reference = Channel.value(tuple([id: "reference"], params.reference ? file(reference) : []))
    ch_jasmine_sample_fai       = Channel.value(tuple([id: "fai"], params.reference ? file("${reference}.fai") : []))
    ch_jasmine_sample_chr_norm  = Channel.value([]) // No chr norm file

    JASMINESV_SAMPLE(
        jasminesv_sample_input,
        ch_jasmine_sample_reference,
        ch_jasmine_sample_fai,
        ch_jasmine_sample_chr_norm
    )

    jasminesv_sample_sources = jasminesv_sample_input
        .map { meta, vcf_list, bams, sample_dists -> tuple(meta.sample ?: meta.id, vcf_list) }

    jasminesv_sample_out_keyed = JASMINESV_SAMPLE.out.vcf
        .map { meta, vcf -> tuple(meta.sample ?: meta.id, tuple(meta, vcf)) }

    jasminesv_sample_out_keyed
        .join(jasminesv_sample_sources)
        .map { sample, leftVal, src_vcfs ->
            def (meta, vcf) = leftVal
            tuple(meta, vcf, src_vcfs)
        } | JASMINE_HEADER_FIX

    ch_jasmine_sample_vcfs = JASMINE_HEADER_FIX.out.vcf
    .map { meta, vcf -> tuple(meta, vcf) }

    sample_filtered = ch_jasmine_sample_vcfs | FILTER_CHR

    sample_filtered | SORT_VCF
    sample_sorted = SORT_VCF.out.vcf

    // ──────────────────────────────────────────────────────────────────────
    // Filter SVs supported by ≥2 callers
    // ──────────────────────────────────────────────────────────────────────

    bcftools_sample_input = sample_sorted
    .map { meta, vcf ->
        def v = vcf.toString()
        def updated_meta = [id: meta.sample ?: meta.id, sample: meta.sample ?: meta.id, step: "caller_support"]
        def idx = file(v + '.csi')
        if( !idx.exists() ) idx = file(v + '.tbi')
        def idx_out = idx.exists() ? idx : []
        tuple(updated_meta, file(v), idx_out)
    }

    CALLER_SUPPORT_FILTER(
        bcftools_sample_input,
        Channel.value([]), // regions
        Channel.value([]), // targets
        Channel.value([])  // samples
    )

    // Consensus summary - simple approach
    consensus_summary_input = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()
        .map { vcf_list -> tuple([id: "consensus_summary"], vcf_list) }

    SUMMARIZE_CALLER_MERGED(
        consensus_summary_input,
        Channel.value("consensus")
    )

    PLOT_CONSENSUS(
        SUMMARIZE_CALLER_MERGED.out.json
            .map { meta, json -> tuple([id: "consensus_plot"], [json]) },
        Channel.value("Consensus SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // SAMPLE LEVEL AF ANNOTATION + FILTERING + ANNOTSV ANNOTATION
    // ──────────────────────────────────────────────────────────────────────

    ch_per_sample_input = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> tuple(meta, vcf) }

    ch_svdb_in_occ  = Channel.value(params.svdb_in_occ ?: [])
    ch_svdb_in_frq  = Channel.value(params.svdb_in_frq ?: [])
    ch_svdb_out_occ = Channel.value(params.svdb_out_occ ?: [])
    ch_svdb_out_frq = Channel.value(params.svdb_out_frq ?: [])
    ch_svdb_dbs     = Channel.value(params.svdb_databases ? params.svdb_databases.collect { file(it) } : [])
    ch_svdb_bedpe   = Channel.value([])

    SVDB_QUERY_SAMPLE(
        ch_per_sample_input,
        ch_svdb_in_occ,
        ch_svdb_in_frq,
        ch_svdb_out_occ,
        ch_svdb_out_frq,
        ch_svdb_dbs,
        ch_svdb_bedpe
    )

    ch_per_sample_bcftools_input = SVDB_QUERY_SAMPLE.out.vcf
    .map { meta, annotated_vcf ->
        def updated_meta = [id: meta.sample ?: meta.id, sample: meta.sample ?: meta.id, step: "af_filter"]
        def idx = file(annotated_vcf.toString() + '.csi')
        if( !idx.exists() ) idx = file(annotated_vcf.toString() + '.tbi')
        def idx_out = idx.exists() ? idx : []
        tuple(updated_meta, file(annotated_vcf), idx_out)
    }

    ch_bcftools_regions = Channel.value([])
    ch_bcftools_targets = Channel.value([])
    ch_bcftools_samples = Channel.value([])

    AF_FILTER(
        ch_per_sample_bcftools_input,
        ch_bcftools_regions,
        ch_bcftools_targets,
        ch_bcftools_samples
    )

    // Filtered summary - simple approach
    filtered_summary_input = AF_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()
        .map { vcf_list -> tuple([id: "filtered_summary"], vcf_list) }

    SUMMARIZE_CALLER_MERGED_FILTERED(
        filtered_summary_input,
        Channel.value("filtered")
    )

    PLOT_FILTERED(
        SUMMARIZE_CALLER_MERGED_FILTERED.out.json
            .map { meta, json -> tuple([id: "filtered_plot"], [json]) },
        Channel.value("Filtered SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // prepare channel for AnnotSV annotations
    // ──────────────────────────────────────────────────────────────────────

    if(!annotsv_annotations) {
        ANNOTSV_INSTALLANNOTATIONS()
        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { [[id:"annotsv"], it] }
            .collect()
            .set { ch_annotsv_annotations }
    } else {
        ch_annotsv_annotations_input = Channel.fromPath(annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        if(annotsv_annotations.endsWith(".tar.gz")){
            UNTAR_ANNOTSV(ch_annotsv_annotations_input)
            UNTAR_ANNOTSV.out.untar
                .collect()
                .set { ch_annotsv_annotations }
        } else {
            ch_annotsv_annotations = Channel.fromPath(annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        }
    }

    ch_candidate_genes      = Channel.value(tuple([id: "candidate_genes"], []))
    ch_false_positive_snv   = Channel.value(tuple([id: "false_positive_snv"], []))
    ch_gene_transcripts     = Channel.value(tuple([id: "gene_transcripts"], []))

    ANNOTSV_PER_SAMPLE_RAW(
        SVDB_QUERY_SAMPLE.out.vcf
            .map { meta, vcf -> 
                def updated_meta = meta.clone()
                updated_meta.id = "${meta.sample ?: meta.id}_annotated"
                tuple(updated_meta, vcf, [], []) 
            },
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    ANNOTSV_PER_SAMPLE(
        AF_FILTER.out.vcf
            .map { meta, vcf ->
                def updated_meta = [
                    id: "${meta.sample ?: meta.id}_filtered",
                    sample: meta.sample ?: meta.id, 
                    step: "final_annotation"
                ]
                tuple(updated_meta, vcf, [], [])
            },
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    // ──────────────────────────────────────────────────────────────────────
    // Continue to cohort-level analyses
    // ──────────────────────────────────────────────────────────────────────

    sample_consensus_vcfs = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()

    jasminesv_cohort_input = sample_consensus_vcfs
        .map { vcf_list ->
        tuple([id: "cohort"], vcf_list, [], [])
    }

    ch_jasmine_cohort_reference = Channel.value(tuple([id: "reference"], params.reference ? file(reference) : []))
    ch_jasmine_cohort_fai       = Channel.value(tuple([id: "fai"], params.reference ? file("${reference}.fai") : []))
    ch_jasmine_cohort_chr_norm  = Channel.value([]) // No chr norm file

    JASMINESV_COHORT(
        jasminesv_cohort_input,
        ch_jasmine_cohort_reference,
        ch_jasmine_cohort_fai,
        ch_jasmine_cohort_chr_norm
    )

    // ──────────────────────────────────────────────────────────────────────
    // SV annotation using SVDB (cohort-level)
    // ──────────────────────────────────────────────────────────────────────

    svdb_cohort_input = JASMINESV_COHORT.out.vcf
        .map { meta, cohort_vcf ->
        tuple(meta, cohort_vcf)
    }

    ch_svdb_cohort_in_occ  = Channel.value(params.svdb_in_occ ?: [])
    ch_svdb_cohort_in_frq  = Channel.value(params.svdb_in_frq ?: [])
    ch_svdb_cohort_out_occ = Channel.value(params.svdb_out_occ ?: [])
    ch_svdb_cohort_out_frq = Channel.value(params.svdb_out_frq ?: [])
    ch_svdb_cohort_dbs     = Channel.value(params.svdb_databases ? params.svdb_databases.collect { file(it) } : [])
    ch_svdb_cohort_bedpe   = Channel.value([])

    SVDB_QUERY_COHORT(
        svdb_cohort_input,
        ch_svdb_cohort_in_occ,
        ch_svdb_cohort_in_frq,
        ch_svdb_cohort_out_occ,
        ch_svdb_cohort_out_frq,
        ch_svdb_cohort_dbs,
        ch_svdb_cohort_bedpe
    )

    ch_candidate_genes_cohort    = Channel.value(tuple([id: "candidate_genes"], []))
    ch_false_positive_snv_cohort = Channel.value(tuple([id: "false_positive_snv"], []))
    ch_gene_transcripts_cohort   = Channel.value(tuple([id: "gene_transcripts"], []))

    ANNOTSV_COHORT_RAW(
        SVDB_QUERY_COHORT.out.vcf
            .map { meta, vcf -> 
                def updated_meta = [
                    id: "${meta.id}_annotated",
                    sample: "${meta.id}_annotated",
                ]
                tuple(updated_meta, vcf, [], [])
            },
        ch_annotsv_annotations,
        ch_candidate_genes_cohort,
        ch_false_positive_snv_cohort,
        ch_gene_transcripts_cohort
    )

    // Cohort summaries - one per cohort directory
    SUMMARIZE_COHORT_ANNOTATED(
        SVDB_QUERY_COHORT.out.vcf
            .map { meta, vcf -> tuple([id: "cohort_annotated_summary"], vcf) },
        Channel.value("cohort_annotated")
    )

    PLOT_COHORT_ANNOTATED(
        SUMMARIZE_COHORT_ANNOTATED.out.json
            .map { meta, json -> tuple([id: "cohort_annotated_plot"], [json]) },
        Channel.value("Cohort Annotated SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // Filter annotated SVs based on AF
    // ──────────────────────────────────────────────────────────────────────

    bcftools_cohort_input = SVDB_QUERY_COHORT.out.vcf.map { meta, annotated_vcf ->
        tuple(meta, annotated_vcf, [])
    }

    AF_FILTER_COHORT(
        bcftools_cohort_input,
        ch_bcftools_regions,
        ch_bcftools_targets,
        ch_bcftools_samples
    )

    // Cohort summaries - one per cohort directory
    SUMMARIZE_COHORT_FILTERED(
        AF_FILTER_COHORT.out.vcf
            .map { meta, vcf -> tuple([id: "cohort_filtered_summary"], vcf) },
        Channel.value("cohort_filtered")
    )

    PLOT_COHORT_FILTERED(
        SUMMARIZE_COHORT_FILTERED.out.json
            .map { meta, json -> tuple([id: "cohort_filtered_plot"], [json]) },
        Channel.value("Cohort Filtered SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // Structural variant annotation using AnnotSV
    // ──────────────────────────────────────────────────────────────────────

    annotsv_input = AF_FILTER_COHORT.out.vcf.map { meta, filtered_vcf ->
        def updated_meta = [
            id: 'cohort_filtered',
            sample: 'cohort_filtered'
        ]
        tuple(updated_meta, filtered_vcf, [], [])
    }

    ANNOTSV_COHORT(
        annotsv_input,
        ch_annotsv_annotations,
        ch_candidate_genes_cohort,
        ch_false_positive_snv_cohort,
        ch_gene_transcripts_cohort
)

    // ──────────────────────────────────────────────────────────────────────
    // Collate and save software versions
    // ──────────────────────────────────────────────────────────────────────

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(CUTESV.out.versions)
    ch_versions = ch_versions.mix(SEVERUS_WITH_CONTROL.out.versions)
    ch_versions = ch_versions.mix(SEVERUS_NO_CONTROL.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_COHORT.out.versions)
    ch_versions = ch_versions.mix(SORT_VCF.out.versions)
    ch_versions = ch_versions.mix(SVDB_QUERY_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(SVDB_QUERY_COHORT.out.versions)
    ch_versions = ch_versions.mix(CALLER_SUPPORT_FILTER.out.versions)
    ch_versions = ch_versions.mix(AF_FILTER.out.versions)
    ch_versions = ch_versions.mix(AF_FILTER_COHORT.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_PER_SAMPLE_RAW.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_PER_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_COHORT.out.versions)

    // Handle conditional AnnotSV installation versions
    if(!annotsv_annotations) {
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'ontvar_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
        multiqc_report         = MULTIQC.out.report.toList()
        versions               = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
