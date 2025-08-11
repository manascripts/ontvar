process SUMMARIZE_SV_COUNTS {
    tag "${meta.id}"

    label 'process_low'
    
    publishDir "${params.outdir}/stats", mode: 'copy', pattern: '*_sv_summary.tsv'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_sv_summary.tsv")

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    echo -e "sample_id\tsv_type\tsv_count" > ${meta.id}_sv_summary.tsv
    if [[ "${vcf}" == *.gz ]]; then
        zcat "${vcf}" | grep -v '^#' | cut -f8 | grep -o 'SVTYPE=\\w\\+' | cut -d= -f2 | sort | uniq -c |
        awk '{print "'${meta.id}'\t" \$2 "\t" \$1}' >> ${meta.id}_sv_summary.tsv
    else
        grep -v '^#' "${vcf}" | cut -f8 | grep -o 'SVTYPE=\\w\\+' | cut -d= -f2 | sort | uniq -c |
        awk '{print "'${meta.id}'\t" \$2 "\t" \$1}' >> ${meta.id}_sv_summary.tsv
    fi
    """
}
