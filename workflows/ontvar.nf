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
include { SNIFFLES } from '../modules/nf-core/sniffles/main'
include { CUTESV   } from '../modules/nf-core/cutesv/main'
include { SEVERUS  } from '../modules/nf-core/severus/main'
include { RENAME_VCF } from '../modules/local/rename_vcf/main'
include { JASMINESV as JASMINESV_SAMPLE } from '../modules/nf-core/jasminesv/main'
include { JASMINESV as JASMINESV_COHORT } from '../modules/nf-core/jasminesv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV    } from '../modules/nf-core/annotsv/annotsv/main'           
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLERS          } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_SAMPLE_MERGED    } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_MERGED    } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_ANNOTATED } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_FILTERED  } from '../modules/local/summarize_sv_counts/main'
include { SVDB_QUERY } from '../modules/nf-core/svdb/query/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SAMPLE } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_COHORT } from '../modules/nf-core/bcftools/view/main'


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

    main:
        ch_multiqc_files = Channel.empty()
        ch_versions = Channel.empty()

    // ──────────────────────────────────────────────────────────────────────
    // Variant calling: match cases with controls and run SV callers
    // ──────────────────────────────────────────────────────────────────────

    //println "Running ONTVAR workflow with samplesheet: ${ch_samplesheet}"
    ch_sample_info = ch_samplesheet
    // println "Running ONTVAR workflow with samplesheet: ${ch_samplesheet}"
    // Channel.fromPath(ch_samplesheet, checkIfExists: true)
    //     .view { ch_samplesheet -> 
    //         println "Samplesheet path: ${ch_samplesheet}" 
    //     }
    // println "Samplesheet loaded: ${ch_samplesheet}"
    // ch_sample_info = file(ch_samplesheet)
    //     .splitCsv(header: true, sep: ',', strip: true)
    //     .view()

    cases = ch_sample_info.filter { it[2] == 'case' }
    controls = ch_sample_info.filter { it[2] == 'control' }

    controls_by_match = controls.map { c -> tuple(c[3], c[1]) }

    sv_input = cases
        .map { c -> tuple(c[3], c[0], c[1]) }  // map to [match_id, sample_id, bam_path]
        .join(controls_by_match, remainder: true)  // join on match_id
        .map { it -> 
            def sample_id = it[1] 
            def case_bam = it[2]
            def control_bam = it[3]  // control_bam comes from join
            tuple(sample_id, case_bam, control_bam ?: null)
        }

    // ──────────────────────────────────────────────────────────────────────
    // SV Calling
    // ──────────────────────────────────────────────────────────────────────

    // SNIFFLES
    sniffles_input = sv_input.map { it ->
        tuple([id: "${it[0]}_sniffles"], it[1], file("${it[1]}.bai"))
    }

    SNIFFLES(
        sniffles_input,                                                         // Input 1: [meta, bam, bai]
        Channel.value(tuple([id: "reference"], file(reference))),               // Input 2: [meta, fasta]
        Channel.value(tuple([id: "tandem"], file(params.sniffles.tandem_repeats))), // Input 3: [meta, tandem_file]
        Channel.value(true),                                                    // Input 4: vcf_output
        Channel.value(false)                                                    // Input 5: snf_output
        
    )
    
    
    // CUTESV
    cutesv_input = sv_input.map { it ->
        tuple([id: "${it[0]}_cutesv"], it[1], file("${it[1]}.bai"))
    }
    
    CUTESV(
        cutesv_input,
        Channel.value(tuple([id: "reference"], file(reference)))
    )

    // SEVERUS  
    severus_input = sv_input.map { it ->
        tuple([id: it[0]], it[1], file("${it[1]}.bai"), it[2] ?: [], it[2] ? file("${it[2]}.bai") : [], [])
    }
    
    SEVERUS(
        severus_input,                                                          // Input 1: [meta, target_bam, target_bai, control_bam, control_bai, vcf]
        Channel.value(tuple([id: "vntr"], file(params.severus.vntr_bed)))       // Input 2: [meta, vntr_bed]
    )

    severus_vcfs  = SEVERUS.out.somatic_vcf.map { meta, vcf -> tuple(meta, vcf, 'severus') } | RENAME_VCF
    
    // Merge all caller outputs and summarize together
    all_caller_vcfs = SNIFFLES.out.vcf.mix(CUTESV.out.vcf, severus_vcfs)
        .map { meta, vcf -> tuple(meta, vcf) }
    all_caller_vcfs | SUMMARIZE_CALLERS

    // ──────────────────────────────────────────────────────────────────────
    // Gather SV caller outputs per sample
    // ──────────────────────────────────────────────────────────────────────

    sv_calls_by_sample = all_caller_vcfs.groupTuple(by: 0)

    // ──────────────────────────────────────────────────────────────────────
    // Run Jasmine to merge SVs from callers per sample  
    // ──────────────────────────────────────────────────────────────────────

    jasminesv_sample_input = sv_calls_by_sample.map { meta, vcf_list -> 
        tuple([id: "${meta}_jasmine", sample: meta], vcf_list, [], [])
    }

    JASMINESV_SAMPLE(
        jasminesv_sample_input,                                                 // Input 1: [meta, vcfs, bams, sample_dists]
        Channel.value(tuple([id: "reference"], params.jasminesv_sample.fasta ? file(reference) : [])),     // Input 2: [meta, fasta]
        Channel.value(tuple([id: "fai"], params.jasminesv_sample.fasta_fai ? file("${reference}.fai") : [])), // Input 3: [meta, fasta_fai]
        Channel.value(params.jasminesv_sample.chr_norm ? file(params.jasminesv_sample.chr_norm) : [])      // Input 4: chr_norm
    )

    // ──────────────────────────────────────────────────────────────────────
    // Filter SVs supported by ≥2 callers
    // ──────────────────────────────────────────────────────────────────────

    bcftools_sample_input = JASMINESV_SAMPLE.out.vcf.map { meta, merged_vcf ->
        // def meta_map = meta instanceof Map && meta.containsKey('id') ? meta : [id: meta]
        def tbi = file("${merged_vcf}.tbi")
        def csi = file("${merged_vcf}.csi")
        def index = tbi.exists() ? tbi : (csi.exists() ? csi : null)
        tuple(meta, merged_vcf, index)
    }

    BCFTOOLS_VIEW_SAMPLE(
        bcftools_sample_input,                                                  // Input 1: [meta, vcf, index]
        Channel.value(params.bcftools_view_sample.regions ? file(params.bcftools_view_sample.regions) : []), // Input 2: regions
        Channel.value(params.bcftools_view_sample.targets ? file(params.bcftools_view_sample.targets) : []), // Input 3: targets
        Channel.value(params.bcftools_view_sample.samples ? file(params.bcftools_view_sample.samples) : [])  // Input 4: samples
    )

    // ──────────────────────────────────────────────────────────────────────
    // SV summary after Jasmine merge
    // ──────────────────────────────────────────────────────────────────────

    BCFTOOLS_VIEW_SAMPLE.out.vcf | SUMMARIZE_SAMPLE_MERGED

    // ──────────────────────────────────────────────────────────────────────
    // Merge sample consensus VCFs into cohort consensus
    // ──────────────────────────────────────────────────────────────────────

    sample_consensus_vcfs = BCFTOOLS_VIEW_SAMPLE.out.vcf.map { meta, vcf -> vcf }.collect()

    jasminesv_cohort_input = sample_consensus_vcfs.map { vcf_list -> 
        tuple([id: "cohort"], vcf_list, [], [])  // meta, vcfs, bams (empty), sample_dists (empty)
    }

    JASMINESV_COHORT(
        jasminesv_cohort_input,                                                 // Input 1: [meta, vcfs, bams, sample_dists]
        Channel.value(tuple([id: "reference"], params.reference ? file(reference) : [])),              // Input 2: [meta, fasta]
        Channel.value(tuple([id: "fai"], params.jasminesv_cohort.fasta_fai ? file("${reference}.fai") : [])), // Input 3: [meta, fasta_fai]
        Channel.value(params.jasminesv_cohort.chr_norm ? file(params.jasminesv_cohort.chr_norm) : [])  // Input 4: chr_norm
    )

    // ──────────────────────────────────────────────────────────────────────
    // SV annotation using SVDB (cohort-level)
    // ──────────────────────────────────────────────────────────────────────

    svdb_input = JASMINESV_COHORT.out.vcf.map { meta, cohort_vcf ->
        // def meta_map = meta instanceof Map && meta.containsKey('id') ? meta : [id: meta]
        tuple(meta, cohort_vcf)
    }

    SVDB_QUERY(
        svdb_input,                                                             // Input 1: [meta, vcf]
        Channel.value(params.svdb_query.in_occ ?: []),                         // Input 2: in_occ
        Channel.value(params.svdb_query.in_frq ?: []),                         // Input 3: in_frq
        Channel.value(params.svdb_query.out_occ ?: []),                        // Input 4: out_occ
        Channel.value(params.svdb_query.out_frq ?: []),                        // Input 5: out_frq
        Channel.value(params.svdb_query.databases ? params.svdb_query.databases.collect { file(it) } : []), // Input 6: databases
        Channel.value([])                                                       // Input 7: bedpe
    )

    // SV summary on annotated cohort VCF
    SVDB_QUERY.out.vcf | SUMMARIZE_COHORT_ANNOTATED

    // ──────────────────────────────────────────────────────────────────────
    // Filter annotated SVs based on AF
    // ──────────────────────────────────────────────────────────────────────

    bcftools_cohort_input = SVDB_QUERY.out.vcf.map { meta, annotated_vcf ->
        // def meta_map = meta instanceof Map && meta.containsKey('id') ? meta : [id: meta]
        tuple(meta, annotated_vcf, [])  // meta, vcf, index (empty)
    }

    BCFTOOLS_VIEW_COHORT(
        bcftools_cohort_input,                                                  // Input 1: [meta, vcf, index]
        Channel.value(params.bcftools_view_cohort.regions ? file(params.bcftools_view_cohort.regions) : []), // Input 2: regions
        Channel.value(params.bcftools_view_cohort.targets ? file(params.bcftools_view_cohort.targets) : []), // Input 3: targets
        Channel.value(params.bcftools_view_cohort.samples ? file(params.bcftools_view_cohort.samples) : [])  // Input 4: samples
    )

    BCFTOOLS_VIEW_COHORT.out.vcf | SUMMARIZE_COHORT_FILTERED

    // ──────────────────────────────────────────────────────────────────────
    // Structural variant annotation using AnnotSV
    // ──────────────────────────────────────────────────────────────────────

    annotsv_input = BCFTOOLS_VIEW_COHORT.out.vcf.map { meta, filtered_vcf ->
        // def meta_map = meta instanceof Map && meta.containsKey('id') ? meta : [id: meta]
        tuple(meta, filtered_vcf, [], [])  // meta, sv_vcf, sv_vcf_index (empty), candidate_small_variants (empty)
    }

    ANNOTSV(
        annotsv_input,                                                          // Input 1: [meta, sv_vcf, sv_vcf_index, candidate_small_variants]
        Channel.value(tuple([id: "annotations"], params.annotsv.annotations ? file(params.annotsv.annotations) : [])),         // Input 2: [meta, annotations]
        Channel.value(tuple([id: "candidate_genes"], params.annotsv.candidate_genes ? file(params.annotsv.candidate_genes) : [])), // Input 3: [meta, candidate_genes]
        Channel.value(tuple([id: "false_positive_snv"], params.annotsv.false_positive_snv ? file(params.annotsv.false_positive_snv) : [])), // Input 4: [meta, false_positive_snv]
        Channel.value(tuple([id: "gene_transcripts"], params.annotsv.gene_transcripts ? file(params.annotsv.gene_transcripts) : []))    // Input 5: [meta, gene_transcripts]
    )

    // ──────────────────────────────────────────────────────────────────────
    // Collate and save software versions
    // ──────────────────────────────────────────────────────────────────────

    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(CUTESV.out.versions)
    ch_versions = ch_versions.mix(SEVERUS.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_COHORT.out.versions)
    ch_versions = ch_versions.mix(SVDB_QUERY.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_COHORT.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV.out.versions)

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_ontvar_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // ──────────────────────────────────────────────────────────────────────
    // MultiQC
    // ──────────────────────────────────────────────────────────────────────

    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),                                             // Input 1: files
        ch_multiqc_config.toList(),                                             // Input 2: config
        ch_multiqc_custom_config.toList(),                                      // Input 3: custom_config
        ch_multiqc_logo.toList(),                                               // Input 4: logo
        [],                                                                     // Input 5: extra_files
        []                                                                      // Input 6: extra_dirs
    )

    emit:
        multiqc_report         = MULTIQC.out.report.toList()                    // channel: /path/to/multiqc_report.html
        versions               = ch_versions                                     // channel: [ path(versions.yml) ]
        cohort_filtered_vcf    = BCFTOOLS_VIEW_COHORT.out.vcf.map { meta, vcf -> vcf }
        cohort_annotated_table = ANNOTSV.out.tsv.map { meta, tsv -> tsv }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
