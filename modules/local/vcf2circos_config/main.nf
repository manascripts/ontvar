process VCF2CIRCOS_CONFIG {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c3ba2dce2296a1b15cbc39792107842b8c89af08aac761395d8861fffef4f9d6/data' :
        'community.wave.seqera.io/library/wget:1.25.0--b191b6c4876b14bb' }"
    

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("*.tar.gz")           , emit: archive
    path "versions.yml"                         , emit: versions

    script:
    """
    wget -q -O config.tar.gz "${url}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | awk '{print \$1":"\$3}')
    END_VERSIONS
    """
}