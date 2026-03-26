process VCF2CIRCOS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/67904c22ad13d9f21fed06139880b4610ed9e15b79d193a1eb533fe0e6f45b2e/data' :
        'community.wave.seqera.io/library/bcftools_vcf2circos:cc324af5624f6e11' }"

    // local container:
    // 05-analysis-tmp/container/vcf2circos1.1_bcftools.sif
    
    input:
    tuple val(meta), path(vcf)
    val assembly
    path config_dir
    path patch

    output:
    tuple val(meta), path("*.*"), emit: plot
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def format = task.ext.format ?: 'html'

    """

    # 1. Use Python to merge your 'patch' and the 'Static' path into one file.
    # This avoids bind mounts and sed escaping issues entirely.
    python3 -c "
import json
# Patch file (default_params.json)
with open('${patch}', 'r') as f:
    params = json.load(f)

# Options.json (with static path)
with open('${config_dir}/Static/options.json', 'r') as f:
    options = json.load(f)

# Merge them: Ensure 'Static' points to the right place 
# and combine with your patched parameters
options['Static'] = '${config_dir}/Static'
merged = {**params, **options}

with open('final_config.json', 'w') as f:
    json.dump(merged, f, indent=4)
"
    # Standardize VCF for the tool
    bcftools norm -m -any ${vcf} -Ov -o split.vcf

    vcf2circos \\
        --input split.vcf \\
        --output ${prefix}.${format} \\
        --assembly ${assembly} \\
        --options final_config.json \\
        ${args}

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		vcf2circos: \$(vcf2circos --version 2>&1 | grep "Version:" | sed 's/Version: //')
		bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
	END_VERSIONS
    """
}