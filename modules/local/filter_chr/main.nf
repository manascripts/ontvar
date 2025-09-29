process FILTER_CHR {
    tag "${meta.sample}"

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${vcf.baseName}.filtered.vcf.gz")

    when:
      task.ext.when == null || task.ext.when

    script:
    """
    if [[ "${vcf}" == *.gz ]]; then
        zcat ${vcf} | awk 'BEGIN{OFS="\\t"} /^#/ {print} /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\\t/ {print}' | bgzip > ${vcf.baseName}.filtered.vcf.gz
    else
        awk 'BEGIN{OFS="\\t"} /^#/ {print} /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\\t/ {print}' ${vcf} | bgzip > ${vcf.baseName}.filtered.vcf.gz
    fi
    """
}
