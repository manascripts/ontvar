process ANNOTSV_TSV_TO_VCF {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/177dd97d4ce9701108c67c61e3d7a5eb522ee80f6c6fc5b262b6da7c61f16d1f/data':
        'community.wave.seqera.io/library/python:3.12.10--22c5b6d5afb90e4f' }"

    input:
    tuple val(meta), path(original_vcf), path(annotsv_tsv)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"                   , emit: versions

    script:
    // Respecting your ID structure: Prefix stays as meta.id
    prefix = task.ext.prefix ?: "${meta.id}"
    def annotation_mode = 'full'

    """
    #!/usr/bin/env python3
    import csv
    import gzip
    import platform
    import os
    import re

    # --- 1. Read AnnotSV TSV ---
    annot_mode = "${annotation_mode}"
    
    # We only skip columns that are redundant to the VCF structure itself
    SKIP_COLS = {
        'AnnotSV_ID', 'SV_chrom', 'SV_start', 'SV_end', 'SV_length', 'SV_type',
        'Samples_ID', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT',
        'Annotation_mode'
    }

    tsv_path = "${annotsv_tsv}"
    annot_by_key = {}
    annot_cols = []

    if os.path.exists(tsv_path):
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            all_cols = reader.fieldnames if reader.fieldnames else []
            
            # CLINICAL FIX: We keep all columns except the VCF structural ones.
            # This ensures Gene_name, ACMG_class, etc., are preserved.
            annot_cols = [c for c in all_cols if c not in SKIP_COLS]

            for row in reader:
                # Ensure we only process 'full' mode lines as requested
                if row.get('Annotation_mode') != annot_mode:
                    continue
                
                chrom = row.get('SV_chrom', '')
                if chrom and not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                
                key = (chrom, row.get('SV_start', ''), row.get('ID', ''))
                if key not in annot_by_key:
                    annot_by_key[key] = row

    # --- 2. Build INFO header lines ---
    def make_vcf_safe(s):
        # Replace characters that break VCF/Circos formats
        return s.replace('(', '_').replace(')', '').replace('/', '_').replace("'", '')\\
                .replace(' ', '_').replace('-', '_').replace('.', '_')

    info_headers = []
    for col in annot_cols:
        safe_id = make_vcf_safe(col)
        info_headers.append(f'##INFO=<ID={safe_id},Number=.,Type=String,Description="AnnotSV: {col}">')

    # --- 3. Process VCF ---
    vcf_path = "${original_vcf}"
    out_path = "${prefix}.vcf"
    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f, open(out_path, 'w') as out_f:
        for line in f:
            # Header handling
            if line.startswith('##'):
                # vcf2circos compatibility: SVLEN and END must be Number=1
                if 'ID=SVLEN,' in line or 'ID=END,' in line:
                    line = re.sub(r'Number=[^,]+', 'Number=1', line)
                out_f.write(line)
                continue
                
            if line.startswith('#CHROM'):
                for h in info_headers:
                    out_f.write(h + '\\n')
                out_f.write(line)
                continue

            # Record handling
            fields = line.split('\\t')
            if len(fields) < 8:
                out_f.write(line)
                continue

            chrom, pos, vid = fields[0], fields[1], fields[2]
            key = (chrom, pos, vid)
            
            if key in annot_by_key:
                row = annot_by_key[key]
                additions = []
                for col in annot_cols:
                    val = row.get(col, '')
                    if val and val not in ['', 'NA', '.', 'nan', 'None']:
                        safe_id = make_vcf_safe(col)
                        # Standard VCF encoding for special characters in values
                        safe_val = str(val).replace(';', '%3B').replace('=', '%3D').replace(',', '%2C').replace(' ', '%20')
                        additions.append(f'{safe_id}={safe_val}')
                
                if additions:
                    fields[7] = fields[7] + ';' + ';'.join(additions)
            
            out_f.write('\\t'.join(fields))

    # --- 4. Versions ---
    with open("versions.yml", 'w') as f:
        f.write(f'"${task.process}":\\n')
        f.write(f'    python: {platform.python_version()}\\n')
    """
}