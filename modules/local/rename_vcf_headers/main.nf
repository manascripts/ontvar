process RENAME_VCF_HEADERS {
    tag "${meta.sample}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    publishDir "${params.outdir}/renamed_vcfs", mode: 'copy', pattern: '*.vcf*'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf)

    when:
    task.ext.when == null || task.ext.when

    script:
     """
    if [[ "${vcf}" == *.gz ]]; then
        zcat ${vcf} | awk 'BEGIN {OFS="\\t"} /^#CHROM/ { \$NF="${meta.sample}" } {print}' | bgzip > tmp.vcf.gz && mv tmp.vcf.gz ${vcf}
    else
        awk 'BEGIN {OFS="\\t"} /^#CHROM/ { \$NF="${meta.sample}" } {print}' ${vcf} > tmp.vcf && mv tmp.vcf ${vcf}
    fi
    """
}
