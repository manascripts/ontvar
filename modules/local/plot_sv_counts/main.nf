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
    from collections import defaultdict

    # Set aesthetic parameters
    plt.rcParams["font.family"] = "Helvetica"
    plt.rcParams["font.size"] = 13.5
    sns.set_style("whitegrid")

    def parse_json_structure(content, source_name):
        data = []
        
        # Handle empty array
        if isinstance(content, list) and len(content) == 0:
            print(f"Warning: {source_name} contains empty array")
            return data
        
        # Get analysis type
        analysis_type = content.get('analysis_type', 'unknown')
        
        if analysis_type == 'multi_sample':
            # Parse multi-sample structure
            samples = content.get('samples', {})
            for sample_id, sample_data in samples.items():
                sample_analysis_type = sample_data.get('analysis_type', 'unknown')
                
                if sample_analysis_type == 'multi_caller':
                    # Extract from combined_stats
                    combined = sample_data.get('combined_stats', {})
                    sv_types = combined.get('sv_types', {})
                    for svtype, type_data in sv_types.items():
                        if isinstance(type_data, dict) and 'count' in type_data:
                            data.append({
                                'source': sample_id,
                                'caller': 'combined',
                                'svtype': svtype,
                                'count': type_data['count']
                            })
                    
                    # Also extract individual caller data
                    callers = sample_data.get('callers', {})
                    for caller_name, caller_data in callers.items():
                        sv_types = caller_data.get('sv_types', {})
                        for svtype, type_data in sv_types.items():
                            if isinstance(type_data, dict) and 'count' in type_data:
                                data.append({
                                    'source': sample_id,
                                    'caller': caller_name,
                                    'svtype': svtype,
                                    'count': type_data['count']
                                })
                
                elif sample_analysis_type == 'single_vcf':
                    # Extract from single VCF structure
                    sv_types = sample_data.get('sv_types', {})
                    for svtype, type_data in sv_types.items():
                        if isinstance(type_data, dict) and 'count' in type_data:
                            data.append({
                                'source': sample_id,
                                'caller': source_name,
                                'svtype': svtype,
                                'count': type_data['count']
                            })
        
        elif analysis_type == 'multi_caller':
            # Parse multi-caller structure
            # Extract from combined_stats
            combined = content.get('combined_stats', {})
            sv_types = combined.get('sv_types', {})
            for svtype, type_data in sv_types.items():
                if isinstance(type_data, dict) and 'count' in type_data:
                    data.append({
                        'source': source_name,
                        'caller': 'combined',
                        'svtype': svtype,
                        'count': type_data['count']
                    })
            
            # Also extract individual caller data
            callers = content.get('callers', {})
            for caller_name, caller_data in callers.items():
                sv_types = caller_data.get('sv_types', {})
                for svtype, type_data in sv_types.items():
                    if isinstance(type_data, dict) and 'count' in type_data:
                        data.append({
                            'source': source_name,
                            'caller': caller_name,
                            'svtype': svtype,
                            'count': type_data['count']
                        })
        
        elif analysis_type == 'single_vcf':
            # Parse single VCF structure
            sv_types = content.get('sv_types', {})
            for svtype, type_data in sv_types.items():
                if isinstance(type_data, dict) and 'count' in type_data:
                    data.append({
                        'source': source_name,
                        'caller': source_name,
                        'svtype': svtype,
                        'count': type_data['count']
                    })
        
        return data

    def load_json_data(json_files):
        all_data = []
        for json_file in json_files:
            source_name = Path(json_file).stem.replace('_summary', '').replace('_sv_counts', '')
            print(f"Processing: {json_file}")
            
            with open(json_file, 'r') as f:
                content = json.load(f)
                parsed_data = parse_json_structure(content, source_name)
                all_data.extend(parsed_data)
        
        return all_data

    def create_empty_plot(prefix, title, message="No variants found"):
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.text(0.5, 0.5, message, 
               ha='center', va='center', fontsize=20, fontweight='bold')
        ax.text(0.5, 0.4, f'Analysis: {title}', 
               ha='center', va='center', fontsize=14, style='italic')
        ax.axis('off')
        plt.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{prefix}_sv_counts.png', dpi=300, bbox_inches='tight')
        plt.close()

    def create_summary_plots(df, prefix, title):
        if df.empty:
            print(f"No data to plot for {title}")
            create_empty_plot(prefix, title)
            return
        
        # Figure 1: Bar plots - Total and by Type
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Total counts by source (using combined caller data)
        if 'source' in df.columns and 'count' in df.columns:
            # Filter to only 'combined' caller to avoid double counting
            df_combined = df[df['caller'] == 'combined'] if 'combined' in df['caller'].values else df
            
            if not df_combined.empty:
                total_counts = df_combined.groupby('source')['count'].sum().reset_index()
                total_counts = total_counts.sort_values('count', ascending=False)
                
                colors = sns.color_palette("Set2", len(total_counts))
                axes[0].bar(range(len(total_counts)), total_counts['count'], color=colors)
                axes[0].set_xticks(range(len(total_counts)))
                axes[0].set_xticklabels(total_counts['source'], rotation=45, ha='right')
                axes[0].set_title('Total SV Counts by Sample', fontweight='bold')
                axes[0].set_ylabel('Count')
                axes[0].grid(axis='y', alpha=0.3)
            else:
                axes[0].text(0.5, 0.5, 'No source data', ha='center', va='center')
                axes[0].axis('off')
        else:
            axes[0].text(0.5, 0.5, 'No source data', ha='center', va='center')
            axes[0].axis('off')
        
        # Plot 2: Counts by SV type (across all sources)
        if 'svtype' in df.columns and 'count' in df.columns:
            df_combined = df[df['caller'] == 'combined'] if 'combined' in df['caller'].values else df
            
            if not df_combined.empty:
                type_counts = df_combined.groupby('svtype')['count'].sum().reset_index()
                type_counts = type_counts.sort_values('count', ascending=False)
                
                colors = sns.color_palette("Set3", len(type_counts))
                axes[1].bar(range(len(type_counts)), type_counts['count'], color=colors)
                axes[1].set_xticks(range(len(type_counts)))
                axes[1].set_xticklabels(type_counts['svtype'], rotation=45, ha='right')
                axes[1].set_title('SV Counts by Type', fontweight='bold')
                axes[1].set_ylabel('Count')
                axes[1].grid(axis='y', alpha=0.3)
            else:
                axes[1].text(0.5, 0.5, 'No type data', ha='center', va='center')
                axes[1].axis('off')
        else:
            axes[1].text(0.5, 0.5, 'No type data', ha='center', va='center')
            axes[1].axis('off')
        
        plt.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{prefix}_sv_counts.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Figure 2: Stacked bar plot with log10 scale (if multiple sources)
        df_combined = df[df['caller'] == 'combined'] if 'combined' in df['caller'].values else df
        
        if 'source' in df_combined.columns and 'svtype' in df_combined.columns:
            if df_combined['source'].nunique() > 1:
                pivot_df = df_combined.pivot_table(values='count', index='source', columns='svtype', fill_value=0)
                
                # Apply log10 transformation (add 1 to avoid log(0))
                pivot_df_log = np.log10(pivot_df + 1)
                
                fig, ax = plt.subplots(figsize=(12, 6))
                pivot_df_log.plot(kind='bar', stacked=True, ax=ax, colormap='Set3')
                
                ax.set_title(f'{title} - Distribution by Type', fontsize=16, fontweight='bold')
                ax.set_xlabel('')  # Remove x-axis title
                ax.set_ylabel('Count (log10 scale)')
                ax.legend(title='SV Type', loc='best')
                ax.grid(axis='y', alpha=0.3)
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                
                plt.savefig(f'{prefix}_sv_counts_stacked.png', dpi=300, bbox_inches='tight')
                plt.close()
        
        # Figure 3: Caller comparison (if multiple callers per sample)
        if 'caller' in df.columns and df['caller'].nunique() > 1:
            # Exclude 'combined' for this plot
            df_callers = df[df['caller'] != 'combined']
            
            if not df_callers.empty and df_callers['caller'].nunique() > 1:
                caller_counts = df_callers.groupby(['source', 'caller'])['count'].sum().reset_index()
                
                # Create grouped bar chart
                fig, ax = plt.subplots(figsize=(14, 6))
                
                sources = caller_counts['source'].unique()
                callers = caller_counts['caller'].unique()
                x = range(len(sources))
                width = 0.8 / len(callers)
                
                colors = sns.color_palette("Set2", len(callers))
                
                for i, caller in enumerate(callers):
                    caller_data = caller_counts[caller_counts['caller'] == caller]
                    counts = [caller_data[caller_data['source'] == s]['count'].sum() 
                             for s in sources]
                    offset = width * (i - len(callers)/2 + 0.5)
                    ax.bar([xi + offset for xi in x], counts, width, 
                          label=caller, color=colors[i])
                
                ax.set_xlabel('')  # Remove x-axis title
                ax.set_ylabel('Count')
                ax.set_title(f'{title} - Caller Comparison', fontsize=16, fontweight='bold')
                ax.set_xticks(x)
                ax.set_xticklabels(sources, rotation=45, ha='right')
                ax.legend(title='Caller', loc='best')
                ax.grid(axis='y', alpha=0.3)
                plt.tight_layout()
                
                plt.savefig(f'{prefix}_sv_counts_callers.png', dpi=300, bbox_inches='tight')
                plt.close()

    # Main execution
    json_files = "${json_files}".split()
    
    # Load data
    data = load_json_data(json_files)
    df = pd.DataFrame(data)
    
    # Debug: print DataFrame info
    print(f"\\nTotal records loaded: {len(df)}")
    if not df.empty:
        print(f"DataFrame columns: {df.columns.tolist()}")
        print(f"DataFrame shape: {df.shape}")
        print(f"\\nSample data:")
        print(df.head(10))
    
    # Create plots
    create_summary_plots(df, '${prefix}', '${plot_title}')
    
    # Create versions file
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        
        import matplotlib
        f.write(f'    matplotlib: {matplotlib.__version__}\\n')
        
        import seaborn
        f.write(f'    seaborn: {seaborn.__version__}\\n')
        
        import pandas
        f.write(f'    pandas: {pandas.__version__}\\n')
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