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
include { ANNOTSV_TSV_TO_VCF as ANNOTSV_TSV_TO_VCF_PER_SAMPLE_RAW } from '../modules/local/annotsv_tsv_to_vcf/main'
include { ANNOTSV_TSV_TO_VCF as ANNOTSV_TSV_TO_VCF_PER_SAMPLE     } from '../modules/local/annotsv_tsv_to_vcf/main'
include { ANNOTSV_TSV_TO_VCF as ANNOTSV_TSV_TO_VCF_COHORT_RAW     } from '../modules/local/annotsv_tsv_to_vcf/main'
include { ANNOTSV_TSV_TO_VCF as ANNOTSV_TSV_TO_VCF_COHORT         } from '../modules/local/annotsv_tsv_to_vcf/main'
include { VCF2CIRCOS_CONFIG } from '../modules/local/vcf2circos_config/main'
include { UNTAR as UNTAR_VCF2CIRCOS_CONFIG } from '../modules/nf-core/untar/main'
include { VCF2CIRCOS as VCF2CIRCOS_SAMPLE      } from '../modules/local/vcf2circos/main'
include { VCF2CIRCOS as VCF2CIRCOS_SAMPLE_RAW  } from '../modules/local/vcf2circos/main'
include { VCF2CIRCOS as VCF2CIRCOS_COHORT_RAW  } from '../modules/local/vcf2circos/main'
include { VCF2CIRCOS as VCF2CIRCOS_COHORT      } from '../modules/local/vcf2circos/main'
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
        vcf2circos_config

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ──────────────────────────────────────────────────────────────────────
    // FASTQ CONCATENATION (if multiple files per sample)
    // ──────────────────────────────────────────────────────────────────────
    
    ch_raw_inputs = ch_samplesheet.map { group, sample, type, input_type, path ->
        def meta = [
            id: group,          // Current process ID (mutable)
            group_id: group,    // Protected anchor (immutable)
            sample: sample, 
            type: type, 
            input_type: input_type
        ]
        return [meta, path]
    }

    // Concatenate FASTQs from directories
    ch_fastq_prep = ch_raw_inputs.filter { it[0].input_type == 'fastq' }
        .map { meta, path ->
            def p = file(path)
            def files = p.isDirectory() ? p.listFiles().findAll { it.name =~ /.(fastq|fq)(\.gz)?$/ }.sort() : [p]
            return [ meta, files ]
        }

    CAT_FASTQ(ch_fastq_prep.map { meta, files -> [ meta + [id: meta.sample], files ] })
    
    // ──────────────────────────────────────────────────────────────────────
    // ALIGNMENT (minimap2 for long-read sequencing)
    // ──────────────────────────────────────────────────────────────────────

    MINIMAP2_ALIGN (
        CAT_FASTQ.out.reads,
        Channel.value([ [id:'ref'], file(params.reference) ]),
        true, 'bai', false, true
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    ch_aligned_bams = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index) 
        .map { meta, bam, bai -> [ meta, bam, bai ] }

    ch_existing_bams = ch_raw_inputs.filter { it[0].input_type == 'bam' }
        .map { meta, path ->
            def bam = file(path)
            def pathStr = path.toString()
            def bai = file("${pathStr}.bai").exists() ? file("${pathStr}.bai") : 
                      file(pathStr.replaceAll(/\.bam$/, ".bai")).exists() ? file(pathStr.replaceAll(/\.bam$/, ".bai")) : null
            if (!bai) { error "Missing BAI for ${meta.id} at ${path}" }
            return [ meta, bam, bai ]
        }

    ch_all_bams = ch_aligned_bams.mix(ch_existing_bams)

    // Map cases and controls with sample_id, then join
    ch_cases    = ch_all_bams.filter { it[0].type == 'case' }
    ch_controls = ch_all_bams.filter { it[0].type == 'control' }

    // Join on Group ID (meta.id)
    ch_sv_input = ch_all_bams
        .map { meta, bam, bai -> [ meta.group_id, meta, bam, bai ] }
        .groupTuple(by: 0) 
        .map { gid, metas, bams, bais ->
            def c_idx = metas.findIndexOf { it.type == 'case' }
            def n_idx = metas.findIndexOf { it.type == 'control' }

            if (c_idx == -1) return null

            def meta = metas[c_idx]
            def cbam = bams[c_idx]
            def cbai = bais[c_idx]
            
            def tbam = (n_idx != -1) ? bams[n_idx] : []
            def tbai = (n_idx != -1) ? bais[n_idx] : []

            return [ meta, cbam, cbai, tbam, tbai ]
        }
        .filter { it != null }

    // ──────────────────────────────────────────────────────────────────────
    // SV Calling
    // ──────────────────────────────────────────────────────────────────────

    // cbam = Case BAM (Tumor) tbam = Control BAM (Normal)

    SNIFFLES (
        ch_sv_input.map { meta, cbam, cbai, tbam, tbai -> [meta, cbam, cbai] },
        [[id:'ref'], file(params.reference)],
        [[id:'tr'], params.tandem_repeats ? file(params.tandem_repeats) : []],
        true, false
    )

    CUTESV (
        ch_sv_input.map { meta, cbam, cbai, tbam, tbai -> [meta, cbam, cbai] },
        [[id:'ref'], file(params.reference)]
    )

    // SEVERUS Routing
    // Channel branch based on whether tbam is an empty list or a file
    ch_sv_input.branch {
        paired: it[3] instanceof Path || it[3] instanceof String
        single: it[3] == []
    }.set { ch_severus_in }

    SEVERUS_WITH_CONTROL (
        ch_severus_in.paired.map { meta, cbam, cbai, tbam, tbai -> 
            [ meta + [id: "${meta.group_id}_tumor_normal"], cbam, cbai, tbam, tbai, [] ] 
        },
        [[id:'vntr'], file(params.vntr_bed)]
    )

    SEVERUS_NO_CONTROL (
        ch_severus_in.single.map { meta, cbam, cbai, tbam, tbai -> 
            [ meta + [id: "${meta.group_id}_tumor_only"], cbam, cbai, [], [], [] ] 
        },
        [[id:'vntr'], file(params.vntr_bed)]
    )

    // ──────────────────────────────────────────────────────────────────────
    // FIX SAMPLE NAMES in VCF HEADERS
    // ──────────────────────────────────────────────────────────────────────

    ch_sniffles_vcf = SNIFFLES.out.vcf | RENAME_VCF_HEADERS_SNIFFLES
    ch_cutesv_vcf   = CUTESV.out.vcf   | RENAME_VCF_HEADERS_CUTESV
    ch_severus_vcf = SEVERUS_WITH_CONTROL.out.somatic_vcf
        .mix(SEVERUS_NO_CONTROL.out.somatic_vcf) 
        | RENAME_VCF_HEADERS_SEVERUS

    // Raw caller summaries
    ch_all_caller_vcfs = ch_sniffles_vcf.mix(ch_cutesv_vcf, ch_severus_vcf)
        .map { meta, vcf -> 
            return [ [id: meta.group_id, group_id: meta.group_id], vcf ] 
        }

    ch_raw_summary_input = ch_sniffles_vcf.mix(ch_cutesv_vcf, ch_severus_vcf)
        .map { meta, vcf -> "${meta.group_id}|${vcf}" } 
        .collect()
        .map { all_pairs -> 
            tuple([id: "raw_calls_summary"], all_pairs) 
        }

    ch_raw_files = ch_all_caller_vcfs
        .map { meta, vcf -> vcf }
        .collect()

    SUMMARIZE_CALLERS(
        ch_raw_summary_input, 
        ch_raw_files,          
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

    sv_calls_by_sample = ch_all_caller_vcfs.groupTuple(by: 0)

    // ──────────────────────────────────────────────────────────────────────
    // Run Jasmine to merge SVs from callers per sample
    // ──────────────────────────────────────────────────────────────────────

    // Prepare Jasmine input channels (per-sample)
    ch_jasmine_sample_reference = Channel.value(tuple([id: "reference"], params.reference ? file(reference) : []))
    ch_jasmine_sample_fai       = Channel.value(tuple([id: "fai"], params.reference ? file("${reference}.fai") : []))
    ch_jasmine_sample_chr_norm  = Channel.value([]) // No chr norm file

    JASMINESV_SAMPLE(
        sv_calls_by_sample.map { meta, vcfs -> 
            [ meta + [id: "${meta.id}_consensus", step: "consensus"], vcfs, [], [] ] 
        },
        ch_jasmine_sample_reference,
        ch_jasmine_sample_fai,
        ch_jasmine_sample_chr_norm
    )
    // .view { meta, list, _1, _2 -> "DEBUG: Sample ${meta.id} has ${list.size()} VCFs" }

    JASMINESV_SAMPLE.out.vcf
        .map { meta, vcf -> [ meta.group_id, meta, vcf ] }
        .join(sv_calls_by_sample.map { meta, vcfs -> [ meta.group_id, vcfs ] }, by: 0)
        .map { gid, meta, vcf, src_vcfs ->
            // This gives JASMINE_HEADER_FIX exactly what it needs:
            // [ meta, merged_vcf, [list_of_original_vcfs] ]
            return [ meta, vcf, src_vcfs ]
        } | JASMINE_HEADER_FIX

        JASMINE_HEADER_FIX.out.vcf | FILTER_CHR | SORT_VCF
        sample_sorted = SORT_VCF.out.vcf

    // ──────────────────────────────────────────────────────────────────────
    // Filter SVs supported by ≥2 callers
    // ──────────────────────────────────────────────────────────────────────

    bcftools_sample_input = sample_sorted
        .map { meta, vcf ->
            def v = vcf.toString()
            def idx = file(v + '.csi').exists() ? file(v + '.csi') : (file(v + '.tbi').exists() ? file(v + '.tbi') : [])
            // Robust meta update
            def updated_meta = meta + [id: "${meta.group_id}_caller_support", step: "caller_support"]
            return [ updated_meta, file(v), idx ]
        }

    CALLER_SUPPORT_FILTER(bcftools_sample_input, [], [], [])
    
    // Consensus summary
    consensus_summary_input = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> "${meta.group_id}|${vcf}" }
        .collect()
        .map { vcf_list -> tuple([id: "consensus_summary"], vcf_list) }

    ch_consensus_files = CALLER_SUPPORT_FILTER.out.vcf.map { meta, vcf -> vcf }.collect()

    SUMMARIZE_CALLER_MERGED(
        consensus_summary_input,
        ch_consensus_files,
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

    // AF Filter following SVDB
    AF_FILTER(
        SVDB_QUERY_SAMPLE.out.vcf.map { meta, vcf ->
            def v = vcf.toString()
            def idx = file(v + '.csi').exists() ? file(v + '.csi') : (file(v + '.tbi').exists() ? file(v + '.tbi') : [])
            return [ meta + [id: "${meta.group_id}_af_filter", step: "af_filter"], vcf, idx ]
        },
        [], [], []
    )

    // Filtered summary - simple approach
    filtered_summary_input = AF_FILTER.out.vcf
        .map { meta, vcf -> "${meta.group_id}|${vcf}" }
        .collect()
        .map { vcf_list -> tuple([id: "filtered_summary"], vcf_list) }

    ch_filtered_files = AF_FILTER.out.vcf.map { meta, vcf -> vcf }.collect()

    SUMMARIZE_CALLER_MERGED_FILTERED(
        filtered_summary_input,
        ch_filtered_files,
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
        SVDB_QUERY_SAMPLE.out.vcf.map { meta, vcf -> 
        [ meta + [id: "${meta.group_id}_raw_annotated"], vcf, [], [] ] 
    },
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    // Join original VCF with AnnotSV TSV output for per-sample raw
    ch_per_sample_raw_for_merge = SVDB_QUERY_SAMPLE.out.vcf
        .map { meta, vcf -> [ meta.group_id, meta, vcf ] }
        .join(ANNOTSV_PER_SAMPLE_RAW.out.tsv.map { meta, tsv -> [ meta.group_id, tsv ] }, by: 0)
        .map { gid, meta, vcf, tsv -> 
            [ meta + [id: "${gid}_raw_annotated"], vcf, tsv ] 
        }

    ANNOTSV_TSV_TO_VCF_PER_SAMPLE_RAW(ch_per_sample_raw_for_merge)

    ANNOTSV_PER_SAMPLE(
        AF_FILTER.out.vcf.map { meta, vcf ->
        [ meta + [id: "${meta.group_id}_filtered_annotated"], vcf, [], [] ]
    },
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    // Join original VCF with AnnotSV TSV output for per-sample filtered
    ch_per_sample_filtered_for_merge = AF_FILTER.out.vcf
        .map { meta, vcf -> [ meta.group_id, meta, vcf ] }
        .join(ANNOTSV_PER_SAMPLE.out.tsv.map { meta, tsv -> [ meta.group_id, tsv ] }, by: 0)
        .map { gid, meta, vcf, tsv -> 
            [ meta + [id: "${gid}_filtered_annotated"], vcf, tsv ] 
        }

    ANNOTSV_TSV_TO_VCF_PER_SAMPLE(ch_per_sample_filtered_for_merge)

    // ──────────────────────────────────────────────────────────────────────
    // Continue to cohort-level analyses
    // ──────────────────────────────────────────────────────────────────────

    jasminesv_cohort_input = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()
        .map { vcf_list -> 
            tuple([id: "cohort", group_id: "cohort", step: "cohort_merge"], vcf_list, [], []) 
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
        SVDB_QUERY_COHORT.out.vcf.map { meta, vcf -> 
            [ meta + [id: "cohort_annotated"], vcf, [], [] ] 
        },
        ch_annotsv_annotations,
        ch_candidate_genes_cohort,
        ch_false_positive_snv_cohort,
        ch_gene_transcripts_cohort
    )

    ch_cohort_raw_for_merge = SVDB_QUERY_COHORT.out.vcf
        .map { meta, vcf -> [ meta.group_id, vcf ] }
        .join(ANNOTSV_COHORT_RAW.out.tsv.map { meta, tsv -> [ meta.group_id, tsv ] })
        .map { gid, vcf, tsv ->
            [ [id: "cohort_annotated", group_id: gid], vcf, tsv ]
        }
    
    ANNOTSV_TSV_TO_VCF_COHORT_RAW(ch_cohort_raw_for_merge)

    ch_cohort_annot_files = SVDB_QUERY_COHORT.out.vcf.map { meta, vcf -> vcf }.collect()

    // Cohort summaries - one per cohort directory
    SUMMARIZE_COHORT_ANNOTATED(
        SVDB_QUERY_COHORT.out.vcf
            .map { meta, vcf -> "cohort|${vcf}" }
            .collect()
            .map { vcf_list -> tuple([id: "cohort_annotated_summary"], vcf_list) },
        ch_cohort_annot_files,
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

    bcftools_cohort_input = SVDB_QUERY_COHORT.out.vcf.map { meta, vcf ->
        def v = vcf.toString()
        def idx = file(v + '.csi').exists() ? file(v + '.csi') : (file(v + '.tbi').exists() ? file(v + '.tbi') : [])
        tuple(meta, vcf, idx)
    }

    AF_FILTER_COHORT(
        bcftools_cohort_input,
        [],
        [],
        []
    )

    ch_cohort_filt_files = AF_FILTER_COHORT.out.vcf.map { meta, vcf -> vcf }.collect()

    // Cohort summaries - one per cohort directory
    SUMMARIZE_COHORT_FILTERED(
        AF_FILTER_COHORT.out.vcf
            .map { meta, vcf -> "cohort|${vcf}" }
            .collect()
            .map { vcf_list -> tuple([id: "cohort_filtered_summary"], vcf_list) },
        ch_cohort_filt_files,
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

    ANNOTSV_COHORT(
        AF_FILTER_COHORT.out.vcf.map { meta, vcf ->
            [ meta + [id: "cohort_filtered"], vcf, [], [] ]
        },
        ch_annotsv_annotations,
        ch_candidate_genes_cohort,
        ch_false_positive_snv_cohort,
        ch_gene_transcripts_cohort
    )

    ch_cohort_filtered_for_merge = AF_FILTER_COHORT.out.vcf
        .map { meta, vcf -> [ meta.group_id, vcf ] }
        .join(ANNOTSV_COHORT.out.tsv.map { meta, tsv -> [ meta.group_id, tsv ] })
        .map { gid, vcf, tsv ->
            [ [id: "cohort_filtered", group_id: gid], vcf, tsv ]
        }

    ANNOTSV_TSV_TO_VCF_COHORT(ch_cohort_filtered_for_merge)


    // ──────────────────────────────────────────────────────────────────────
    // CIRCOS PLOTS - now use TSV-TO-VCF output (valid VCF headers)
    // ──────────────────────────────────────────────────────────────────────
    if (params.plot_circos) {
        def circos_sources = (params.circos_sources ?: []) as List
        ch_assembly = Channel.value(params.genome_assembly ?: 'hg38')
        ch_vcf_patch = Channel.value(file("${projectDir}/assets/vcf2circos/default_params.json"))

        // Resolve vcf2circos config
        if (vcf2circos_config) {
            if (vcf2circos_config.toString().endsWith('.tar.gz')) {
                ch_vcf2circos_tarball = Channel.fromPath(vcf2circos_config).map { [[id: 'vcf2circos_config_user'], it] }
                UNTAR_VCF2CIRCOS_CONFIG(ch_vcf2circos_tarball)
                ch_vcf2circos_config_dir = UNTAR_VCF2CIRCOS_CONFIG.out.untar.map { meta, dir -> dir }.collect()
            } else {
                ch_vcf2circos_config_dir = Channel.fromPath(vcf2circos_config, type: 'dir').collect()
            }
        } else {
            VCF2CIRCOS_CONFIG(Channel.value([[id: 'vcf2circos_config'], params.vcf2circos_url]))
            UNTAR_VCF2CIRCOS_CONFIG(VCF2CIRCOS_CONFIG.out.archive)
            ch_vcf2circos_config_dir = UNTAR_VCF2CIRCOS_CONFIG.out.untar.map { meta, dir -> dir }.collect()
        }

        // Pass the resolved directory channel directly
        if (circos_sources.contains('per_sample')) {
            VCF2CIRCOS_SAMPLE_RAW(ANNOTSV_TSV_TO_VCF_PER_SAMPLE_RAW.out.vcf, ch_assembly, ch_vcf2circos_config_dir, ch_vcf_patch)
        }
        if (circos_sources.contains('per_sample_filtered')) {
            VCF2CIRCOS_SAMPLE(ANNOTSV_TSV_TO_VCF_PER_SAMPLE.out.vcf, ch_assembly, ch_vcf2circos_config_dir, ch_vcf_patch)
        }
        if (circos_sources.contains('cohort_annotated')) {
            VCF2CIRCOS_COHORT_RAW(ANNOTSV_TSV_TO_VCF_COHORT_RAW.out.vcf, ch_assembly, ch_vcf2circos_config_dir, ch_vcf_patch)
        }
        if (circos_sources.contains('cohort_filtered')) {
            VCF2CIRCOS_COHORT(ANNOTSV_TSV_TO_VCF_COHORT.out.vcf, ch_assembly, ch_vcf2circos_config_dir, ch_vcf_patch)
        }
    }

    // ──────────────────────────────────────────────────────────────────────
    // Collate and save software versions
    // ──────────────────────────────────────────────────────────────────────

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(CUTESV.out.versions)
    ch_versions = ch_versions.mix(SEVERUS_WITH_CONTROL.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_COHORT.out.versions)
    ch_versions = ch_versions.mix(SORT_VCF.out.versions)
    ch_versions = ch_versions.mix(SVDB_QUERY_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(CALLER_SUPPORT_FILTER.out.versions)
    ch_versions = ch_versions.mix(AF_FILTER_COHORT.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_COHORT.out.versions)

    // Handle conditional AnnotSV installation versions
    if(!annotsv_annotations) {
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)
    }

    if (params.plot_circos) {
        def circos_sources2 = (params.circos_sources ?: []) as List
        if (circos_sources2.contains('cohort_filtered')) {
            ch_versions = ch_versions.mix(VCF2CIRCOS_COHORT.out.versions)
        }
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
