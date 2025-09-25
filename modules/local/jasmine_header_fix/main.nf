
process JASMINE_HEADER_FIX {
    tag "${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02b45b077ea80378d94e41519f4f2fadf8af69a9f1211e1a227fccd2dd2ee9ac/data':
        'community.wave.seqera.io/library/bcftools_coreutils_htslib_awk_bash:8371466e4ddc03d8' }"

    /*
     - input:
       * meta: map with sample metadata
       * jasmine_vcf: path to jasmine VCF (vcf or vcf.gz)
       * src_vcfs: a Groovy list (val) of source caller VCF paths
     - output: emits channel name "vcf" as used in workflow
    */
    input:
        tuple val(meta), path(jasmine_vcf), val(src_vcfs)

    output:
        tuple val(meta), path("${meta.sample}.jasmine.fixed.vcf.gz"), emit: vcf

    when:
      task.ext.when == null || task.ext.when
        
    script:
    """
    set -euo pipefail

    # Force all temp operations to use current working directory
    export TMPDIR=\$PWD
    export TMP=\$PWD
    export TEMP=\$PWD
    export BCFTOOLS_PLUGINS=""
    
    # Create a local temp directory
    mkdir -p \$PWD/tmp
    export TMPDIR=\$PWD/tmp

    in_vcf='${jasmine_vcf}'
    out_vcf='${meta.sample}.jasmine.fixed.vcf.gz'
    hdr_tmp=all_headers.tmp
    jasmine_hdr=jasmine_headers.tmp
    hdr_union=hdr_union.tmp
    header_add=header_add.tmp

    > \${hdr_tmp}

    # collect ## lines from source VCFs (individual callers)
    for f in ${src_vcfs.collect{ "\"${it}\"" }.join(' ')} ; do
        if [[ -s "\${f}" ]]; then
            if [[ "\${f}" == *.gz ]]; then
                zcat "\${f}" | awk '/^##/ {print}' >> \${hdr_tmp} || true
            else
                awk '/^##/ {print}' "\${f}" >> \${hdr_tmp} || true
            fi
        fi
    done

    # extract header from jasmine vcf
    if [[ "\${in_vcf}" == *.gz ]]; then
        zcat "\${in_vcf}" | awk '/^##/ {print}' > \${jasmine_hdr} || true
    else
        awk '/^##/ {print}' "\${in_vcf}" > \${jasmine_hdr} || true
    fi

    # build union of header lines (dedupe exact lines)
    if [[ -s \${hdr_tmp} ]]; then
        sort -u \${hdr_tmp} > \${hdr_union}
    else
        > \${hdr_union}
    fi

    # compute lines present in union but missing in jasmine header
    if [[ -s \${hdr_union} ]]; then
        comm -23 <(sort \${hdr_union}) <(sort \${jasmine_hdr}) > \${header_add} || true
    else
        > \${header_add}
    fi

    # ensure header_add exists
    if [[ ! -s \${header_add} ]]; then
        echo -n > \${header_add}
    fi

    # annotate jasmine vcf with missing header lines AND sort
    # Use explicit temp directory for bcftools
    if [[ -s \${header_add} ]]; then
        bcftools annotate --header-lines \${header_add} "\${in_vcf}" | \\
        bcftools sort -T \$PWD/tmp -Oz -o "\${out_vcf}"
    else
        bcftools sort -T \$PWD/tmp -Oz -o "\${out_vcf}" "\${in_vcf}"
    fi

    # Index the sorted VCF
    tabix -p vcf "\${out_vcf}"
    """
}
