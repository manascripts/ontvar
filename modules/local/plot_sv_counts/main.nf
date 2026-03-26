process PLOT_SV_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e4ea18715a7792c201df25d215712e88c33d4c6583729cd800d41f6f231a2f95/data' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_seaborn:ea777521476506ba' }"

    input:
    tuple val(meta), path(json_files)
    val(plot_title)

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3

    import json
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import sys
    import os

    # Set aesthetic parameters
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 12, 'font.family': 'sans-serif'})

    def parse_new_json(json_file):
        rows = []
        with open(json_file, 'r') as f:
            content = json.load(f)
            
        stage = content.get('stage', 'unknown')
        samples = content.get('samples', {})

        for sample_id, sample_data in samples.items():
            callers = sample_data.get('callers', {})
            # If no callers (empty file), add a placeholder
            if not callers:
                 rows.append({
                    'Sample': sample_id, 'Caller': 'none',
                    'SVType': 'NONE', 'Count': 0, 'Stage': stage
                })
                 continue

            for caller_name, caller_stats in callers.items():
                sv_types = caller_stats.get('sv_types', {})
                if not sv_types:
                    rows.append({
                        'Sample': sample_id, 'Caller': caller_name,
                        'SVType': 'NONE', 'Count': 0, 'Stage': stage
                    })
                for sv_type, type_data in sv_types.items():
                    rows.append({
                        'Sample': sample_id,
                        'Caller': caller_name,
                        'SVType': sv_type,
                        'Count': type_data.get('count', 0),
                        'Stage': stage
                    })
        return rows

    input_files = "${json_files}".replace('[', '').replace(']', '').split(',')
    # Clean up whitespace and handle potential single vs multiple inputs
    json_files = [f.strip() for f in input_files if f.strip()]
    
    all_rows = []
    for f in json_files:
        if os.path.exists(f):
            all_rows.extend(parse_new_json(f))

    df = pd.DataFrame(all_rows)

    if df.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No variants found in input JSONs", ha='center')
        plt.savefig('${prefix}_sv_counts.png')
        # Touch dummy files to satisfy Nextflow outputs
        Path('${prefix}_sv_counts_callers.png').touch()
        Path('${prefix}_sv_counts_stacked.png').touch()
        sys.exit(0)

    # 1. Total SV Counts per Sample (Stacked by Type)
    # We use 'consensus' if available (Merged stage), otherwise sum callers (Raw stage)
    has_consensus = 'consensus' in df['Caller'].values
    plot_df = df[df['Caller'] == 'consensus'] if has_consensus else df
    
    pivot_df = plot_df.groupby(['Sample', 'SVType'])['Count'].sum().unstack(fill_value=0)
    
    fig, ax = plt.subplots(figsize=(12, 7))
    pivot_df.plot(kind='bar', stacked=True, ax=ax, colormap='viridis')
    ax.set_title('${plot_title}')
    ax.set_ylabel('Number of Variants')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('${prefix}_sv_counts.png', dpi=300)

    # 2. Caller Comparison Plot
    if df['Caller'].nunique() > 1:
        plt.figure(figsize=(12, 7))
        caller_comp = df.groupby(['Sample', 'Caller'])['Count'].sum().unstack(fill_value=0)
        caller_comp.plot(kind='bar', ax=plt.gca())
        plt.title('Comparison of Callers per Sample')
        plt.ylabel('Total SVs Detected')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('${prefix}_sv_counts_callers.png', dpi=300)
    else:
        Path('${prefix}_sv_counts_callers.png').touch()

    # 3. Log-Scale Stacked Plot
    plt.figure(figsize=(12, 7))
    log_pivot = np.log10(pivot_df + 1)
    log_pivot.plot(kind='bar', stacked=True, ax=plt.gca(), colormap='plasma')
    plt.title('${plot_title} (Log10 Scale)')
    plt.ylabel('log10(Count + 1)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('${prefix}_sv_counts_stacked.png', dpi=300)

    # Versions file
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
        f.write(f'    matplotlib: {matplotlib.__version__}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sv_counts.png
    touch ${prefix}_sv_counts_stacked.png
    touch ${prefix}_sv_counts_callers.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}