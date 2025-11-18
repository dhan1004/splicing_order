#!/usr/bin/env python3
"""
Analyze RNA secondary structure features and their relationship to splicing order.

Usage:
    python analyze_structure.py structure_features.tsv

Generates:
    - figures/: Directory with all analysis plots
    - correlation_summary.tsv: Table of correlation statistics
    - structure_analysis_report.txt: Text summary of findings
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
import argparse
import os
import sys
from pathlib import Path

# Set plotting style
sns.set_style('whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

def load_data(filepath):
    """Load structure features TSV file."""
    print(f"\nLoading data from {filepath}...")
    df = pd.read_csv(filepath, sep='\t')
    
    print(f"  Loaded {len(df):,} intron pairs")
    print(f"  Samples: {df['sample_id'].nunique()}")
    print(f"  Genes: {df['gene_id'].nunique()}")
    print(f"  Total reads: {df['total_reads'].sum():,}")
    
    return df

def plot_splicing_distribution(df, output_dir):
    """Plot distribution of splicing order."""
    print("\n1. Analyzing splicing order distribution...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    axes[0].hist(df['fraction_downstream'], bins=50, edgecolor='black', alpha=0.7)
    axes[0].axvline(0.5, color='red', linestyle='--', linewidth=2, label='Random (0.5)')
    axes[0].set_xlabel('Fraction Downstream Spliced First')
    axes[0].set_ylabel('Number of Intron Pairs')
    axes[0].set_title('Distribution of Splicing Order')
    axes[0].legend()
    
    # Cumulative distribution
    sorted_vals = np.sort(df['fraction_downstream'])
    axes[1].plot(sorted_vals, np.arange(1, len(sorted_vals)+1)/len(sorted_vals), linewidth=2)
    axes[1].axvline(0.5, color='red', linestyle='--', linewidth=2, label='Random (0.5)')
    axes[1].set_xlabel('Fraction Downstream Spliced First')
    axes[1].set_ylabel('Cumulative Fraction')
    axes[1].set_title('Cumulative Distribution')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/01_splicing_order_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Statistical test
    t_stat, p_val = stats.ttest_1samp(df['fraction_downstream'], 0.5)
    
    return {
        'mean_fraction': df['fraction_downstream'].mean(),
        'median_fraction': df['fraction_downstream'].median(),
        'ttest_t': t_stat,
        'ttest_p': p_val
    }

def plot_mfe_analysis(df, output_dir):
    """Analyze MFE (thermodynamic stability) effects."""
    print("\n2. Analyzing MFE effects...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # MFE distributions
    axes[0, 0].hist(df['intron1_mfe'].dropna(), bins=50, alpha=0.5, label='Intron 1', edgecolor='black')
    axes[0, 0].hist(df['intron2_mfe'].dropna(), bins=50, alpha=0.5, label='Intron 2', edgecolor='black')
    axes[0, 0].set_xlabel('MFE (kcal/mol)')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('MFE Distribution by Intron')
    axes[0, 0].legend()
    
    # MFE difference distribution
    axes[0, 1].hist(df['mfe_difference'].dropna(), bins=50, edgecolor='black', alpha=0.7)
    axes[0, 1].axvline(0, color='red', linestyle='--', linewidth=2)
    axes[0, 1].set_xlabel('MFE Difference (Intron1 - Intron2)')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Distribution of MFE Differences')
    
    # Scatter: MFE difference vs splicing order
    valid_data = df[['mfe_difference', 'fraction_downstream', 'total_reads']].dropna()
    axes[1, 0].scatter(valid_data['mfe_difference'], valid_data['fraction_downstream'], 
                       s=np.log10(valid_data['total_reads'])*20, alpha=0.3, c=valid_data['total_reads'], 
                       cmap='viridis', edgecolors='black', linewidth=0.5)
    axes[1, 0].axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 0].axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 0].set_xlabel('MFE Difference (Intron1 - Intron2)')
    axes[1, 0].set_ylabel('Fraction Downstream Spliced First')
    axes[1, 0].set_title('MFE Difference vs Splicing Order')
    
    # Correlation
    r_pearson, p_pearson = None, None
    if len(valid_data) > 0:
        r_pearson, p_pearson = pearsonr(valid_data['mfe_difference'], valid_data['fraction_downstream'])
        r_spearman, p_spearman = spearmanr(valid_data['mfe_difference'], valid_data['fraction_downstream'])
        axes[1, 0].text(0.05, 0.95, f'Pearson r={r_pearson:.3f}, p={p_pearson:.2e}\\nSpearman ρ={r_spearman:.3f}, p={p_spearman:.2e}',
                       transform=axes[1, 0].transAxes, verticalalignment='top', 
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Binned analysis
    bins = pd.qcut(df['mfe_difference'].dropna(), q=10, duplicates='drop')
    binned_data = df.groupby(bins).agg({
        'fraction_downstream': 'mean',
        'mfe_difference': 'mean'
    }).dropna()
    
    axes[1, 1].scatter(binned_data['mfe_difference'], binned_data['fraction_downstream'], 
                       s=200, alpha=0.7, edgecolors='black', linewidth=2)
    axes[1, 1].axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 1].axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 1].set_xlabel('MFE Difference (Binned Mean)')
    axes[1, 1].set_ylabel('Mean Fraction Downstream')
    axes[1, 1].set_title('Binned Analysis (10 bins)')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/02_mfe_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return {'mfe_r': r_pearson, 'mfe_p': p_pearson}

def plot_accessibility_analysis(df, output_dir):
    """Analyze splice site accessibility effects."""
    print("\n3. Analyzing splice site accessibility...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    results = {}
    
    # 5' splice site accessibility
    valid_5ss = df[['ss5_accessibility_difference', 'fraction_downstream', 'total_reads']].dropna()
    
    axes[0, 0].scatter(valid_5ss['ss5_accessibility_difference'], valid_5ss['fraction_downstream'],
                       s=np.log10(valid_5ss['total_reads'])*20, alpha=0.3, c=valid_5ss['total_reads'],
                       cmap='viridis', edgecolors='black', linewidth=0.5)
    axes[0, 0].axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 0].axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 0].set_xlabel("5'SS Accessibility Difference (Intron1 - Intron2)")
    axes[0, 0].set_ylabel('Fraction Downstream Spliced First')
    axes[0, 0].set_title("5' Splice Site Accessibility vs Splicing Order")
    
    if len(valid_5ss) > 0:
        r, p = pearsonr(valid_5ss['ss5_accessibility_difference'], valid_5ss['fraction_downstream'])
        results['ss5_r'] = r
        results['ss5_p'] = p
        axes[0, 0].text(0.05, 0.95, f'r={r:.3f}, p={p:.2e}',
                       transform=axes[0, 0].transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 3' splice site accessibility
    valid_3ss = df[['ss3_accessibility_difference', 'fraction_downstream', 'total_reads']].dropna()
    
    axes[0, 1].scatter(valid_3ss['ss3_accessibility_difference'], valid_3ss['fraction_downstream'],
                       s=np.log10(valid_3ss['total_reads'])*20, alpha=0.3, c=valid_3ss['total_reads'],
                       cmap='viridis', edgecolors='black', linewidth=0.5)
    axes[0, 1].axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 1].axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 1].set_xlabel("3'SS Accessibility Difference (Intron1 - Intron2)")
    axes[0, 1].set_ylabel('Fraction Downstream Spliced First')
    axes[0, 1].set_title("3' Splice Site Accessibility vs Splicing Order")
    
    if len(valid_3ss) > 0:
        r, p = pearsonr(valid_3ss['ss3_accessibility_difference'], valid_3ss['fraction_downstream'])
        results['ss3_r'] = r
        results['ss3_p'] = p
        axes[0, 1].text(0.05, 0.95, f'r={r:.3f}, p={p:.2e}',
                       transform=axes[0, 1].transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Combined accessibility
    df['avg_accessibility_diff'] = (df['ss5_accessibility_difference'] + df['ss3_accessibility_difference']) / 2
    valid_avg = df[['avg_accessibility_diff', 'fraction_downstream', 'total_reads']].dropna()
    
    axes[0, 2].scatter(valid_avg['avg_accessibility_diff'], valid_avg['fraction_downstream'],
                       s=np.log10(valid_avg['total_reads'])*20, alpha=0.3, c=valid_avg['total_reads'],
                       cmap='viridis', edgecolors='black', linewidth=0.5)
    axes[0, 2].axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 2].axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 2].set_xlabel('Average Accessibility Difference')
    axes[0, 2].set_ylabel('Fraction Downstream Spliced First')
    axes[0, 2].set_title('Average Splice Site Accessibility')
    
    if len(valid_avg) > 0:
        r, p = pearsonr(valid_avg['avg_accessibility_diff'], valid_avg['fraction_downstream'])
        axes[0, 2].text(0.05, 0.95, f'r={r:.3f}, p={p:.2e}',
                       transform=axes[0, 2].transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Binned analyses (bottom row)
    for idx, (col, title) in enumerate([('ss5_accessibility_difference', "5'SS"),
                                          ('ss3_accessibility_difference', "3'SS"),
                                          ('avg_accessibility_diff', 'Average')]):
        bins = pd.qcut(df[col].dropna(), q=10, duplicates='drop')
        binned = df.groupby(bins).agg({
            'fraction_downstream': 'mean',
            col: 'mean'
        }).dropna()
        
        axes[1, idx].scatter(binned[col], binned['fraction_downstream'],
                            s=200, alpha=0.7, edgecolors='black', linewidth=2)
        axes[1, idx].axhline(0.5, color='red', linestyle='--', linewidth=1, alpha=0.5)
        axes[1, idx].axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
        axes[1, idx].set_xlabel(f'{title} Accessibility Difference (Binned)')
        axes[1, idx].set_ylabel('Mean Fraction Downstream')
        axes[1, idx].set_title(f'{title} Binned Analysis')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/03_accessibility_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return results

def plot_categorical_comparisons(df, output_dir):
    """Compare categorical structure features."""
    print("\n4. Analyzing categorical comparisons...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    results = {}
    
    # More structured intron comparison
    if 'intron1_more_structured' in df.columns and df['intron1_more_structured'].notna().sum() > 0:
        grouped = df.groupby('intron1_more_structured')['fraction_downstream'].apply(list)
        
        positions = [0, 1]
        data_to_plot = [grouped.get(0, []), grouped.get(1, [])]
        
        bp = axes[0].boxplot(data_to_plot, positions=positions, widths=0.6,
                             patch_artist=True, showmeans=True,
                             medianprops=dict(color='red', linewidth=2),
                             meanprops=dict(marker='D', markerfacecolor='blue', markersize=8))
        
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        
        axes[0].axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        axes[0].set_xticks([0, 1])
        axes[0].set_xticklabels(['Intron2 More Structured', 'Intron1 More Structured'])
        axes[0].set_ylabel('Fraction Downstream Spliced First')
        axes[0].set_title('Splicing Order by Structure Complexity')
        
        # Statistical test
        if len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
            u_stat, p_val = mannwhitneyu(data_to_plot[0], data_to_plot[1])
            results['structure_p'] = p_val
            axes[0].text(0.5, 0.95, f'Mann-Whitney U test\\np={p_val:.2e}',
                        transform=axes[0].transAxes, horizontalalignment='center',
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # More accessible intron comparison
    if 'intron1_more_accessible' in df.columns and df['intron1_more_accessible'].notna().sum() > 0:
        grouped = df.groupby('intron1_more_accessible')['fraction_downstream'].apply(list)
        
        positions = [0, 1]
        data_to_plot = [grouped.get(0, []), grouped.get(1, [])]
        
        bp = axes[1].boxplot(data_to_plot, positions=positions, widths=0.6,
                             patch_artist=True, showmeans=True,
                             medianprops=dict(color='red', linewidth=2),
                             meanprops=dict(marker='D', markerfacecolor='blue', markersize=8))
        
        for patch in bp['boxes']:
            patch.set_facecolor('lightgreen')
        
        axes[1].axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        axes[1].set_xticks([0, 1])
        axes[1].set_xticklabels(['Intron2 More Accessible', 'Intron1 More Accessible'])
        axes[1].set_ylabel('Fraction Downstream Spliced First')
        axes[1].set_title('Splicing Order by Splice Site Accessibility')
        
        # Statistical test
        if len(data_to_plot[0]) > 0 and len(data_to_plot[1]) > 0:
            u_stat, p_val = mannwhitneyu(data_to_plot[0], data_to_plot[1])
            results['accessibility_p'] = p_val
            axes[1].text(0.5, 0.95, f'Mann-Whitney U test\\np={p_val:.2e}',
                        transform=axes[1].transAxes, horizontalalignment='center',
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/04_categorical_comparisons.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return results

def generate_correlation_table(df, output_file):
    """Generate correlation summary table."""
    print("\n5. Generating correlation summary table...")
    
    summary_results = []
    
    test_vars = [
        ('MFE Difference', 'mfe_difference'),
        ("5'SS Accessibility Diff", 'ss5_accessibility_difference'),
        ("3'SS Accessibility Diff", 'ss3_accessibility_difference'),
        ('Intron1 Length', 'intron1_length'),
        ('Intron2 Length', 'intron2_length'),
        ('Intron1 GC Content', 'intron1_gc_content'),
        ('Intron2 GC Content', 'intron2_gc_content'),
        ('Intron1 Pairing Fraction', 'intron1_pairing_fraction'),
        ('Intron2 Pairing Fraction', 'intron2_pairing_fraction'),
    ]
    
    for name, col in test_vars:
        if col in df.columns:
            valid = df[[col, 'fraction_downstream']].dropna()
            if len(valid) > 0:
                r_p, p_p = pearsonr(valid[col], valid['fraction_downstream'])
                r_s, p_s = spearmanr(valid[col], valid['fraction_downstream'])
                
                summary_results.append({
                    'Feature': name,
                    'N': len(valid),
                    'Pearson_r': r_p,
                    'Pearson_p': p_p,
                    'Spearman_rho': r_s,
                    'Spearman_p': p_s
                })
    
    summary_df = pd.DataFrame(summary_results)
    summary_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"  Saved correlation table to {output_file}")
    return summary_df

def generate_text_report(df, dist_results, output_file):
    """Generate text summary report."""
    print("\n6. Generating text summary report...")
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("RNA SECONDARY STRUCTURE ANALYSIS - SUMMARY REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Dataset overview
        f.write("1. DATASET OVERVIEW\n")
        f.write(f"   Total intron pairs: {len(df):,}\n")
        f.write(f"   Unique samples: {df['sample_id'].nunique()}\n")
        f.write(f"   Unique genes: {df['gene_id'].nunique()}\n")
        f.write(f"   Total reads: {df['total_reads'].sum():,}\n")
        f.write(f"   Mean reads/pair: {df['total_reads'].mean():.1f}\n")
        f.write(f"   Median reads/pair: {df['total_reads'].median():.0f}\n\n")
        
        # Splicing order distribution
        f.write("2. SPLICING ORDER DISTRIBUTION\n")
        f.write(f"   Mean fraction downstream: {dist_results['mean_fraction']:.3f}\n")
        f.write(f"   Median fraction downstream: {dist_results['median_fraction']:.3f}\n")
        f.write(f"   T-test vs 0.5: t={dist_results['ttest_t']:.3f}, p={dist_results['ttest_p']:.2e}\n")
        
        upstream_bias = (df['fraction_downstream'] < 0.4).sum()
        downstream_bias = (df['fraction_downstream'] > 0.6).sum()
        neutral = ((df['fraction_downstream'] >= 0.4) & (df['fraction_downstream'] <= 0.6)).sum()
        
        f.write(f"   Upstream-biased (<0.4): {upstream_bias:,} ({upstream_bias/len(df)*100:.1f}%)\n")
        f.write(f"   Downstream-biased (>0.6): {downstream_bias:,} ({downstream_bias/len(df)*100:.1f}%)\n")
        f.write(f"   Neutral (0.4-0.6): {neutral:,} ({neutral/len(df)*100:.1f}%)\n\n")
        
        # Key correlations
        f.write("3. STRUCTURE-SPLICING CORRELATIONS\n")
        test_vars = [
            ('MFE Difference', 'mfe_difference'),
            ("5'SS Accessibility", 'ss5_accessibility_difference'),
            ("3'SS Accessibility", 'ss3_accessibility_difference'),
        ]
        
        for name, col in test_vars:
            if col in df.columns:
                valid = df[[col, 'fraction_downstream']].dropna()
                if len(valid) > 0:
                    r, p = pearsonr(valid[col], valid['fraction_downstream'])
                    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
                    f.write(f"   {name:25s} r={r:7.3f}, p={p:.2e} {sig}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("Analysis complete. See figures/ directory for visualizations.\n")
        f.write("="*80 + "\n")
    
    print(f"  Saved text report to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze RNA structure features and splicing order'
    )
    parser.add_argument('input_file', help='Structure features TSV file')
    parser.add_argument('--output-dir', default='figures', 
                       help='Output directory for figures (default: figures)')
    
    args = parser.parse_args()
    
    # Check input file
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file not found: {args.input_file}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("\n" + "="*80)
    print("RNA SECONDARY STRUCTURE ANALYSIS")
    print("="*80)
    
    # Load data
    df = load_data(args.input_file)
    
    # Run analyses
    dist_results = plot_splicing_distribution(df, args.output_dir)
    mfe_results = plot_mfe_analysis(df, args.output_dir)
    acc_results = plot_accessibility_analysis(df, args.output_dir)
    cat_results = plot_categorical_comparisons(df, args.output_dir)
    
    # Generate summary outputs
    corr_table = generate_correlation_table(df, 'correlation_summary.tsv')
    generate_text_report(df, dist_results, 'structure_analysis_report.txt')
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  - {args.output_dir}/*.png (analysis figures)")
    print(f"  - correlation_summary.tsv (correlation statistics)")
    print(f"  - structure_analysis_report.txt (summary report)")
    print("\n")

if __name__ == "__main__":
    main()