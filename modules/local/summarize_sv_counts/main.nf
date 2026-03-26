process SUMMARIZE_SV_COUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198b15b0581f19cfa29c5ff506b138aeb1bedb226d2ca308c46705bab133c1f0/data' :
        'community.wave.seqera.io/library/bcftools_coreutils_gawk_gzip_pruned:e1a91ca0c5f22302' }"

    input:
    tuple val(meta), val(vcf_metadata)
    path(vcf_files, stageAs: "inputs_?/*") 
    val stage_name

    output:
    tuple val(meta), path("${stage_name}_summary.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def vcf_metadata_str = vcf_metadata.join(' ')
    """
    python3 << 'EOF'
import sys, json, statistics, os, gzip, glob
from collections import defaultdict, OrderedDict

def open_vcf(vcf_file):
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    return open(vcf_file, 'r')

def analyze_vcf(vcf_file):
    result = {"total_variants": 0, "sv_types": defaultdict(lambda: {"count": 0, "svlen_data": []})}
    if not os.path.exists(vcf_file):
        return result
    try:
        with open_vcf(vcf_file) as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 8: continue
                result["total_variants"] += 1
                info = fields[7]
                info_dict = dict(item.split('=', 1) for item in info.split(';') if '=' in item)
                svtype = info_dict.get('SVTYPE', info_dict.get('TYPE', 'UNK'))
                if svtype == 'UNK' and fields[4].startswith('<'):
                    svtype = fields[4].strip('<>')
                result["sv_types"][svtype]["count"] += 1
                svlen = None
                if 'SVLEN' in info_dict:
                    try: svlen = abs(int(info_dict['SVLEN'].split(',')[0]))
                    except: pass
                elif 'END' in info_dict:
                    try: svlen = abs(int(info_dict['END']) - int(fields[1]))
                    except: pass
                if svlen is not None and svlen > 0:
                    result["sv_types"][svtype]["svlen_data"].append(svlen)
    except Exception as e:
        print(f"Error reading {vcf_file}: {e}")

    for svtype, data in result["sv_types"].items():
        if data["svlen_data"]:
            lengths = data["svlen_data"]
            data["svlen_min"], data["svlen_max"] = min(lengths), max(lengths)
            data["svlen_mean"] = round(statistics.mean(lengths), 2)
            data["svlen_median"] = round(statistics.median(lengths), 1)
        if "svlen_data" in data: del data["svlen_data"]
    return result

# 1. Identify all files staged in the inputs_X/ folders
staged_paths = glob.glob("inputs_*/*")

# 2. Parse the metadata string
raw_metadata = "${vcf_metadata_str}".split()
sample_to_files = defaultdict(list)

for item in raw_metadata:
    if '|' in item:
        sample_id, original_path = item.split('|', 1)
        filename = os.path.basename(original_path)
    else:
        filename = os.path.basename(item)
        sample_id = filename.split('_')[0]

    # Find which staged path matches this filename
    match = None
    for p in staged_paths:
        if os.path.basename(p) == filename:
            match = p
            # Remove from list so we don't assign the same physical file 
            # to multiple metadata entries if names are identical
            staged_paths.remove(p)
            break
    
    if match:
        sample_to_files[sample_id].append(match)

# 3. Assemble JSON
final_output = OrderedDict([
    ("stage", "${stage_name}"),
    ("analysis_type", "cohort" if len(sample_to_files) > 1 else "individual"),
    ("samples", {})
])

for sample_id, vcf_list in sample_to_files.items():
    sample_entry = {"callers": {}, "summary_stats": {"total_variants": 0}}
    for vcf in vcf_list:
        v_low = vcf.lower()
        if 'sniffles' in v_low:   caller = 'sniffles'
        elif 'cutesv' in v_low:  caller = 'cutesv'
        elif 'severus' in v_low: caller = 'severus'
        elif 'jasmine' in v_low or 'consensus' in v_low: caller = 'consensus'
        else: caller = 'unknown'
        
        stats = analyze_vcf(vcf)
        sample_entry["callers"][caller] = stats
        sample_entry["summary_stats"]["total_variants"] += stats["total_variants"]
    final_output["samples"][sample_id] = sample_entry

with open("${stage_name}_summary.json", "w") as f:
    json.dump(final_output, f, indent=2)
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}