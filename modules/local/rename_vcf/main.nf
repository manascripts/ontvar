process RENAME_VCF {
    tag "${meta.id}_${suffix}"
    label 'process_single'

    input:
      tuple val(meta), path(vcf), val(suffix)
    
    output:
      tuple val(meta), path("${meta.id}_${suffix}.vcf")
    
    when:
      task.ext.when == null || task.ext.when

    script:
      """
      cp ${vcf} ${meta.id}_${suffix}.vcf
      """
}