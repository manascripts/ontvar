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
    // Smart detection of input type
    def is_list = vcf_input instanceof List || (vcf_input instanceof String && vcf_input.contains(' '))
    def vcf_files = is_list ? (vcf_input instanceof List ? vcf_input.join(' ') : vcf_input) : vcf_input.toString()
    """
    python3 << 'EOF'
import sys, json, statistics, subprocess, os
from collections import defaultdict, OrderedDict

def analyze_vcf(vcf_file):
    if not os.path.exists(vcf_file) or os.path.getsize(vcf_file) == 0:
        print(f"Warning: VCF file {vcf_file} is empty or does not exist")
        return {"total_variants": 0, "sv_types": {}}
        
    # Check VCF header for available fields
    header_cmd = f"bcftools view -h {vcf_file}"
    try:
        header_result = subprocess.run(header_cmd, shell=True, capture_output=True, text=True, timeout=30)
        if header_result.returncode != 0:
            print(f"Warning: Could not read header from {vcf_file}")
            return {"total_variants": 0, "sv_types": {}}
    except Exception as e:
        print(f"Warning: Error reading {vcf_file}: {e}")
        return {"total_variants": 0, "sv_types": {}}
    
    # Determine available INFO fields
    has_chr2 = '##INFO=<ID=CHR2' in header_result.stdout
    has_end = '##INFO=<ID=END' in header_result.stdout
    has_svlen = '##INFO=<ID=SVLEN' in header_result.stdout
    has_mateid = '##INFO=<ID=MATEID' in header_result.stdout
    has_event = '##INFO=<ID=EVENT' in header_result.stdout
    
    print(f"Analyzing {os.path.basename(vcf_file)}: CHR2={has_chr2}, END={has_end}, SVLEN={has_svlen}, MATEID={has_mateid}, EVENT={has_event}")
    
    # Build query format based on available fields
    query_fields = ['%CHROM', '%POS', '%ID', '%INFO/SVTYPE']
    if has_chr2:
        query_fields.append('%INFO/CHR2')
    if has_end:
        query_fields.append('%INFO/END')
    if has_svlen:
        query_fields.append('%INFO/SVLEN')
    if has_mateid:
        query_fields.append('%INFO/MATEID')
    if has_event:
        query_fields.append('%INFO/EVENT')
    
    query_format = '\\t'.join(query_fields) + '\\n'
    
    # Query VCF data
    query_cmd = f"bcftools query -f '{query_format}' {vcf_file}"
    try:
        query_result = subprocess.run(query_cmd, shell=True, capture_output=True, text=True, timeout=60)
        if query_result.returncode != 0:
            print(f"Warning: Could not query {vcf_file}")
            return {"total_variants": 0, "sv_types": {}}
    except Exception as e:
        print(f"Warning: Error querying {vcf_file}: {e}")
        return {"total_variants": 0, "sv_types": {}}
    
    sv_counts = defaultdict(int)
    svlen_dict = defaultdict(list)
    processed_events = set()  # Track processed breakend pairs
    bnd_pairs = {}           # Track BND mate relationships
    processed_bnd_ids = set() # Track processed BND IDs
    
    total_lines = 0
    deduplicated_count = 0
    
    for line in query_result.stdout.strip().split('\\n'):
        if not line:
            continue
        total_lines += 1
        fields = line.split('\\t')
        if len(fields) < 4:
            continue
            
        chrom = fields[0]
        pos = int(fields[1]) if fields[1].isdigit() else 0
        variant_id = fields[2]
        svtype = fields[3]
        
        if not svtype or svtype == '.':
            continue
        
        # Parse additional fields based on what's available
        field_idx = 4
        chr2 = fields[field_idx] if has_chr2 and field_idx < len(fields) else "."
        field_idx += 1 if has_chr2 else 0
        
        end = fields[field_idx] if has_end and field_idx < len(fields) else "."
        field_idx += 1 if has_end else 0
        
        svlen = fields[field_idx] if has_svlen and field_idx < len(fields) else None
        field_idx += 1 if has_svlen else 0
        
        mateid = fields[field_idx] if has_mateid and field_idx < len(fields) else "."
        field_idx += 1 if has_mateid else 0
        
        event_id = fields[field_idx] if has_event and field_idx < len(fields) else "."
        
        # Enhanced breakend/translocation handling
        if svtype in {"BND", "TRA"}:
            should_count = True
            
            # Strategy 1: Use EVENT ID if available (JASMINE often provides this)
            if event_id != "." and event_id:
                if event_id in processed_events:
                    should_count = False
                    deduplicated_count += 1
                else:
                    processed_events.add(event_id)
                    
            # Strategy 2: Use MATEID for BND pairs
            elif svtype == "BND" and mateid != "." and mateid:
                if variant_id in processed_bnd_ids or mateid in processed_bnd_ids:
                    should_count = False
                    deduplicated_count += 1
                else:
                    processed_bnd_ids.add(variant_id)
                    processed_bnd_ids.add(mateid)
                    svtype = "TRA"  # Count BND pairs as TRA
                    
            # Strategy 3: Use position-based deduplication for TRA
            elif svtype == "TRA" and chr2 != "." and end != ".":
                try:
                    end_pos = int(end)
                    # Create a canonical key (always smaller chromosome/position first)
                    if chrom < chr2 or (chrom == chr2 and pos < end_pos):
                        key = (chrom, pos, chr2, end_pos)
                    else:
                        key = (chr2, end_pos, chrom, pos)
                    
                    if key in processed_events:
                        should_count = False
                        deduplicated_count += 1
                    else:
                        processed_events.add(key)
                except (ValueError, TypeError):
                    # Fallback: count as single event
                    pass
                    
            # Strategy 4: For other cases, use ID-based deduplication
            else:
                # Check if this looks like a mate pair ID (ends with _1, _2, etc.)
                if variant_id.endswith(('_1', '_2', '.1', '.2')):
                    base_id = variant_id[:-2]
                    if base_id in processed_events:
                        should_count = False
                        deduplicated_count += 1
                    else:
                        processed_events.add(base_id)
            
            if should_count:
                sv_counts[svtype] += 1
                
        else:
            # Standard SV types (DEL, DUP, INS, INV)
            sv_counts[svtype] += 1
        
        # Collect SVLEN data (with filtering for unrealistic values)
        if svlen not in (None, ".", ""):
            try:
                svlen_val = abs(int(float(svlen)))
                # Filter out unrealistic SVLEN values for BND/TRA
                if svtype in {"BND", "TRA"}:
                    # For translocations, very large SVLEN values are often artifacts
                    # Only include if SVLEN seems reasonable (< 1MB for most cases)
                    if svlen_val < 1000000:  # 1MB cutoff
                        svlen_dict[svtype].append(svlen_val)
                else:
                    # For other SV types, include all positive SVLEN values
                    if svlen_val > 0:
                        svlen_dict[svtype].append(svlen_val)
            except (ValueError, TypeError):
                pass
    
    # Build result structure
    sv_types = {}
    for svtype in sorted(sv_counts.keys()):  # Sort for consistent ordering
        entry = {"count": sv_counts[svtype]}
        svlens = svlen_dict.get(svtype, [])
        if svlens:
            entry["svlen_min"] = min(svlens)
            entry["svlen_max"] = max(svlens)
            entry["svlen_mean"] = round(statistics.mean(svlens), 2)
            entry["svlen_median"] = statistics.median(svlens)
            entry["svlen_stdev"] = round(statistics.stdev(svlens) if len(svlens) > 1 else 0, 2)
        else:
            entry.update({
                "svlen_min": None, "svlen_max": None, "svlen_mean": None,
                "svlen_median": None, "svlen_stdev": None
            })
        sv_types[svtype] = entry
    
    total_variants = sum(sv_counts.values())
    print(f"Processed {os.path.basename(vcf_file)}: {total_lines} raw variants -> {total_variants} counted variants")
    if deduplicated_count > 0:
        print(f"  - Deduplicated {deduplicated_count} breakend/translocation events")
    
    return {
        "total_variants": total_variants,
        "sv_types": sv_types
    }

# Smart processing based on input type
vcf_files_str = "${vcf_files}"
stage = "${stage_name}"

# Determine if this is multi-caller or single VCF analysis
vcf_list = [f.strip() for f in vcf_files_str.split() if f.strip()]
is_multi_caller = len(vcf_list) > 1

print(f"Processing stage: {stage}")
print(f"Analysis mode: {'multi_caller' if is_multi_caller else 'single_vcf'}")
print(f"VCF files to process: {len(vcf_list)}")

if is_multi_caller:
    # Multi-caller analysis (for raw calls)
    result = OrderedDict([
        ("stage", stage),
        ("analysis_type", "multi_caller"),
        ("callers", {}),
        ("combined_stats", {
            "total_variants": 0,
            "sv_types": defaultdict(lambda: {"count": 0, "svlen_data": []})
        })
    ])
    
    # Analyze each VCF file (caller)
    for vcf_file in vcf_list:
        if not vcf_file or not os.path.exists(vcf_file):
            print(f"Warning: VCF file {vcf_file} not found")
            continue
        
        # Extract caller name from filename
        basename = os.path.basename(vcf_file)
        if 'sniffles' in basename.lower():
            caller = 'sniffles'
        elif 'cutesv' in basename.lower():
            caller = 'cutesv'
        elif 'severus' in basename.lower():
            caller = 'severus'
        elif 'jasmine' in basename.lower() or 'consensus' in basename.lower():
            caller = 'consensus'
        else:
            # Fallback: use filename without extension
            caller = basename.split('.')[0].split('_')[-1]
        
        print(f"\\nProcessing {caller}: {basename}")
        
        # Analyze this VCF with enhanced deduplication
        stats = analyze_vcf(vcf_file)
        result["callers"][caller] = stats
        
        # Add to combined stats
        result["combined_stats"]["total_variants"] += stats["total_variants"]
        for svtype, data in stats["sv_types"].items():
            result["combined_stats"]["sv_types"][svtype]["count"] += data["count"]
            if data.get("svlen_min") is not None:
                # Collect all SVLEN statistics for later aggregation
                svlen_values = []
                if data.get("svlen_min") is not None:
                    svlen_values.append(data["svlen_min"])
                if data.get("svlen_max") is not None:
                    svlen_values.append(data["svlen_max"])
                if data.get("svlen_mean") is not None:
                    svlen_values.append(data["svlen_mean"])
                if data.get("svlen_median") is not None:
                    svlen_values.append(data["svlen_median"])
                result["combined_stats"]["sv_types"][svtype]["svlen_data"].extend(svlen_values)
    
    # Finalize combined stats
    final_combined = {"total_variants": result["combined_stats"]["total_variants"], "sv_types": {}}
    for svtype in sorted(result["combined_stats"]["sv_types"].keys()):
        data = result["combined_stats"]["sv_types"][svtype]
        entry = {"count": data["count"]}
        svlen_data = [x for x in data["svlen_data"] if x is not None and x > 0]
        if svlen_data:
            entry.update({
                "svlen_min": min(svlen_data),
                "svlen_max": max(svlen_data),
                "svlen_mean": round(statistics.mean(svlen_data), 2),
                "svlen_median": statistics.median(svlen_data),
                "svlen_stdev": round(statistics.stdev(svlen_data) if len(svlen_data) > 1 else 0, 2)
            })
        else:
            entry.update({
                "svlen_min": None, "svlen_max": None, "svlen_mean": None,
                "svlen_median": None, "svlen_stdev": None
            })
        final_combined["sv_types"][svtype] = entry
    
    result["combined_stats"] = final_combined

else:
    # Single VCF analysis (for consensus stages) - CONSISTENT ORDERING
    vcf_file = vcf_list[0] if vcf_list else ""
    print(f"\\nProcessing single VCF: {os.path.basename(vcf_file) if vcf_file else 'None'}")
    
    analyzed_data = analyze_vcf(vcf_file)
    
    # Build result with consistent field ordering using OrderedDict
    result = OrderedDict([
        ("stage", stage),
        ("analysis_type", "single_vcf"),
        ("total_variants", analyzed_data["total_variants"]),
        ("sv_types", analyzed_data["sv_types"])
    ])

# Write output with consistent formatting
with open("${stage_name}_summary.json", "w") as f:
    json.dump(result, f, indent=2)
    
print(f"\\n=== SUMMARY ===")
print(f"Generated summary for stage: {stage}")
print(f"Analysis type: {'multi_caller' if is_multi_caller else 'single_vcf'}")
print(f"Total variants: {result.get('total_variants', 0)}")

if is_multi_caller:
    print("Per-caller breakdown:")
    for caller, stats in result.get("callers", {}).items():
        print(f"  {caller}: {stats.get('total_variants', 0)} variants")
else:
    print("SV type breakdown:")
    for svtype, data in result.get("sv_types", {}).items():
        print(f"  {svtype}: {data.get('count', 0)} variants")

print("===============")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}