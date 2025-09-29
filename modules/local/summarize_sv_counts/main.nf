process SUMMARIZE_SV_COUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/198b15b0581f19cfa29c5ff506b138aeb1bedb226d2ca308c46705bab133c1f0/data' :
        'community.wave.seqera.io/library/bcftools_coreutils_gawk_gzip_pruned:e1a91ca0c5f22302' }"

    input:
    tuple val(meta), path(vcf_input)  // Can be single VCF or list of VCFs
    val stage_name

    output:
    tuple val(meta), path("${stage_name}_summary.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def is_list = vcf_input instanceof List
    def vcf_files = is_list ? vcf_input.join(' ') : vcf_input.toString()
    """
    python3 << 'EOF'
import sys, json, statistics, subprocess, os
from collections import defaultdict, OrderedDict

def analyze_vcf(vcf_file):
    result = {
        "total_variants": 0,
        "sv_types": defaultdict(lambda: {"count": 0, "svlen_data": []})
    }

    if not vcf_file or not os.path.exists(vcf_file):
        return result

    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\\t')
                if len(fields) < 8:
                    continue

                result["total_variants"] += 1
                info = fields[7]

                svtype = "UNK"
                for item in info.split(';'):
                    if item.startswith('SVTYPE='):
                        svtype = item.split('=')[1]
                        break

                result["sv_types"][svtype]["count"] += 1

                for item in info.split(';'):
                    if item.startswith('SVLEN='):
                        try:
                            svlen = abs(int(item.split('=')[1]))
                            result["sv_types"][svtype]["svlen_data"].append(svlen)
                        except:
                            pass
                        break

    except Exception as e:
        print(f"Error processing {vcf_file}: {e}")

    # Compute statistics
    for svtype, data in result["sv_types"].items():
        if data["svlen_data"]:
            lengths = data["svlen_data"]
            data["svlen_min"] = min(lengths)
            data["svlen_max"] = max(lengths)
            data["svlen_mean"] = round(statistics.mean(lengths), 2)
            data["svlen_median"] = round(statistics.median(lengths), 1)
            data["svlen_stdev"] = round(statistics.stdev(lengths) if len(lengths) > 1 else 0, 2)
        del data["svlen_data"]

    return result

# Get inputs
vcf_files_str = "${vcf_files}"
stage = "${stage_name}"
vcf_list = [f.strip() for f in vcf_files_str.split() if f.strip()]

# Detect if we should group by sample based on filenames
sample_data = {}
for vcf_file in vcf_list:
    # Extract sample ID from filename
    basename = os.path.basename(vcf_file)
    if '_' in basename:
        sample_id = basename.split('_')[0]  # e.g., SAMPLE_caller.vcf
    else:
        sample_id = basename.split('.')[0]   # fallback to full basename without extension

    if sample_id not in sample_data:
        sample_data[sample_id] = []
    sample_data[sample_id].append(vcf_file)

# Only use sample grouping if we have multiple samples or multiple VCFs per sample
use_sample_grouping = len(sample_data) > 1 or any(len(vcfs) > 1 for vcfs in sample_data.values())

if use_sample_grouping:
    # Multi-sample nested analysis
    result = OrderedDict([
        ("stage", stage),
        ("analysis_type", "multi_sample"),
        ("samples", {})
    ])

    for sample_id, sample_vcfs in sample_data.items():
        print(f"Processing sample: {sample_id}")
        sample_vcfs_list = sample_vcfs if isinstance(sample_vcfs, list) else [sample_vcfs]

        if len(sample_vcfs_list) > 1:
            # Multi-caller for this sample
            sample_result = OrderedDict([
                ("analysis_type", "multi_caller"),
                ("callers", {}),
                ("combined_stats", {
                    "total_variants": 0,
                    "sv_types": defaultdict(lambda: {"count": 0, "svlen_data": []})
                })
            ])

            caller_map = {}
            for vcf_file in sample_vcfs_list:
                vcf_name = os.path.basename(str(vcf_file)).lower()
                if 'sniffles' in vcf_name:
                    caller_map[vcf_file] = 'sniffles'
                elif 'cutesv' in vcf_name:
                    caller_map[vcf_file] = 'cutesv'
                elif 'severus' in vcf_name:
                    caller_map[vcf_file] = 'severus'
                else:
                    caller_map[vcf_file] = 'unknown'

            for vcf_file, caller in caller_map.items():
                analyzed_data = analyze_vcf(str(vcf_file))
                sample_result["callers"][caller] = analyzed_data

                sample_result["combined_stats"]["total_variants"] += analyzed_data["total_variants"]
                for svtype, data in analyzed_data["sv_types"].items():
                    sample_result["combined_stats"]["sv_types"][svtype]["count"] += data["count"]
                    if "svlen_data" in data:
                        sample_result["combined_stats"]["sv_types"][svtype]["svlen_data"].extend(data["svlen_data"])

            # Compute combined stats
            for svtype, data in sample_result["combined_stats"]["sv_types"].items():
                if data["svlen_data"]:
                    lengths = data["svlen_data"]
                    data["svlen_min"] = min(lengths)
                    data["svlen_max"] = max(lengths)
                    data["svlen_mean"] = round(statistics.mean(lengths), 2)
                    data["svlen_median"] = round(statistics.median(lengths), 1)
                    data["svlen_stdev"] = round(statistics.stdev(lengths) if len(lengths) > 1 else 0, 2)
                    del data["svlen_data"]

        else:
            # Single VCF for this sample
            analyzed_data = analyze_vcf(str(sample_vcfs_list[0]))
            sample_result = OrderedDict([
                ("analysis_type", "single_vcf"),
                ("total_variants", analyzed_data["total_variants"]),
                ("sv_types", analyzed_data["sv_types"])
            ])

        result["samples"][sample_id] = sample_result

else:
    # Regular analysis (cohort-level or legacy)
    if len(vcf_list) > 1:
        # Multi-caller analysis
        result = OrderedDict([
            ("stage", stage),
            ("analysis_type", "multi_caller"),
            ("callers", {}),
            ("combined_stats", {
                "total_variants": 0,
                "sv_types": defaultdict(lambda: {"count": 0, "svlen_data": []})
            })
        ])

        caller_map = {}
        for vcf_file in vcf_list:
            vcf_name = os.path.basename(vcf_file).lower()
            if 'sniffles' in vcf_name:
                caller_map[vcf_file] = 'sniffles'
            elif 'cutesv' in vcf_name:
                caller_map[vcf_file] = 'cutesv'
            elif 'severus' in vcf_name:
                caller_map[vcf_file] = 'severus'
            else:
                caller_map[vcf_file] = 'unknown'

        for vcf_file, caller in caller_map.items():
            analyzed_data = analyze_vcf(vcf_file)
            result["callers"][caller] = analyzed_data

            result["combined_stats"]["total_variants"] += analyzed_data["total_variants"]
            for svtype, data in analyzed_data["sv_types"].items():
                result["combined_stats"]["sv_types"][svtype]["count"] += data["count"]
                if "svlen_data" in data:
                    result["combined_stats"]["sv_types"][svtype]["svlen_data"].extend(data["svlen_data"])

        # Compute combined statistics
        for svtype, data in result["combined_stats"]["sv_types"].items():
            if data["svlen_data"]:
                lengths = data["svlen_data"]
                data["svlen_min"] = min(lengths)
                data["svlen_max"] = max(lengths)
                data["svlen_mean"] = round(statistics.mean(lengths), 2)
                data["svlen_median"] = round(statistics.median(lengths), 1)
                data["svlen_stdev"] = round(statistics.stdev(lengths) if len(lengths) > 1 else 0, 2)
                del data["svlen_data"]

    else:
        # Single VCF analysis
        vcf_file = vcf_list[0] if vcf_list else ""
        analyzed_data = analyze_vcf(vcf_file)

        result = OrderedDict([
            ("stage", stage),
            ("analysis_type", "single_vcf"),
            ("total_variants", analyzed_data["total_variants"]),
            ("sv_types", analyzed_data["sv_types"])
        ])

# Write output
with open("${stage_name}_summary.json", "w") as f:
    json.dump(result, f, indent=2)

print(f"Generated summary for stage: {stage}")
if 'samples' in result:
    print(f"Samples processed: {len(result['samples'])}")
    for sample_id, sample_data in result['samples'].items():
        if 'total_variants' in sample_data:
            print(f"  {sample_id}: {sample_data['total_variants']} variants")
        elif 'combined_stats' in sample_data:
            print(f"  {sample_id}: {sample_data['combined_stats']['total_variants']} variants")
else:
    print(f"Total variants: {result.get('total_variants', 0)}")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
