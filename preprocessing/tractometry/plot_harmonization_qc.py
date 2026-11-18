"""
QC visualization for ComBat harmonization approaches.
Compares raw vs. whole-brain vs. per-tract vs. per-node ComBat
on reference (control) subjects across cohorts.

For each tract and scalar, generates:
  - Row 1: Group mean profiles with 95% CI by cohort
  - Row 2: Site-level mean profiles (to check site effect removal)

Columns: Raw | Whole-brain | Per-tract | Per-node
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
import sys

base = os.environ["PROJECT_DIR"] + "/hcp_tractometry"
out_dir = f'{base}/figures/harmonization_qc'
os.makedirs(out_dir, exist_ok=True)

# Load all datasets
print("Loading datasets...")
datasets = {}
files = {
    'Raw': f'{base}/combined_profiles.csv',
    'Whole-brain': f'{base}/combined_profiles_harmonized.csv',
    'Per-tract': f'{base}/combined_profiles_harmonized_pertract.csv',
    'Per-node': f'{base}/combined_profiles_harmonized_pernode.csv',
}

for label, path in files.items():
    if os.path.exists(path):
        df = pd.read_csv(path, low_memory=False)
        # Only keep reference subjects (controls)
        df = df[df['combat_reference'] == 'yes'].copy()
        datasets[label] = df
        print(f"  {label}: {df['subject'].nunique()} subjects")
    else:
        print(f"  {label}: FILE NOT FOUND ({path})")

approach_names = list(datasets.keys())

# Tracts and scalars to plot
tracts_to_plot = sorted(datasets['Raw']['tractID'].unique())
scalars_to_plot = ['dti_fa', 'dti_md']

scalar_labels = {
    'dti_fa': 'Fractional Anisotropy (FA)',
    'dti_md': 'Mean Diffusivity (MD)',
    'dti_rd': 'Radial Diffusivity (RD)',
    'dti_ad': 'Axial Diffusivity (AD)',
}

# Cohort colors
cohort_colors = {
    'HCP-YA': '#ff7f0e',
    'HCP-Aging': '#1f77b4',
    'Clinical-3T': '#2ca02c',
}

# Plot 1: Cohort mean profiles for each tract/scalar across approaches
print("\nGenerating cohort comparison plots...")
for scalar in scalars_to_plot:
    for tract in tracts_to_plot:
        fig, axes = plt.subplots(2, len(approach_names), figsize=(5*len(approach_names), 8),
                                 sharey='row')
        fig.suptitle(f'{tract} — {scalar_labels.get(scalar, scalar)}', fontsize=14, fontweight='bold')

        for col, approach in enumerate(approach_names):
            df = datasets[approach]
            tract_df = df[df['tractID'] == tract].copy()
            tract_df['nodeID'] = pd.to_numeric(tract_df['nodeID'], errors='coerce')
            tract_df[scalar] = pd.to_numeric(tract_df[scalar], errors='coerce')

            # Top row: individual traces (spaghetti) by cohort
            ax_top = axes[0, col]
            for cohort in sorted(tract_df['cohort'].unique()):
                cohort_data = tract_df[tract_df['cohort'] == cohort]
                color = cohort_colors.get(cohort, 'gray')
                # Plot individual traces (subsample for readability)
                subjects = cohort_data['subject'].unique()
                n_plot = min(50, len(subjects))
                rng = np.random.RandomState(42)
                sample_subs = rng.choice(subjects, n_plot, replace=False)
                for sub in sample_subs:
                    sub_data = cohort_data[cohort_data['subject'] == sub].sort_values('nodeID')
                    ax_top.plot(sub_data['nodeID'], sub_data[scalar],
                               color=color, alpha=0.15, linewidth=0.5)

            ax_top.set_title(approach, fontsize=12)
            if col == 0:
                ax_top.set_ylabel(scalar_labels.get(scalar, scalar), fontsize=10)
            ax_top.set_xlabel('Node')

            # Add cohort legend to first column only
            if col == 0:
                for cohort in sorted(cohort_colors.keys()):
                    if cohort in tract_df['cohort'].unique():
                        ax_top.plot([], [], color=cohort_colors[cohort], linewidth=2, label=cohort)
                ax_top.legend(fontsize=8, loc='upper left')

            # Bottom row: group means with 95% CI
            ax_bot = axes[1, col]
            for cohort in sorted(tract_df['cohort'].unique()):
                cohort_data = tract_df[tract_df['cohort'] == cohort]
                color = cohort_colors.get(cohort, 'gray')
                grouped = cohort_data.groupby('nodeID')[scalar]
                mean = grouped.mean()
                sem = grouped.sem()
                nodes = mean.index

                ax_bot.plot(nodes, mean.values, color=color, linewidth=2, label=cohort)
                ax_bot.fill_between(nodes, (mean - 1.96*sem).values, (mean + 1.96*sem).values,
                                    color=color, alpha=0.2)

            if col == 0:
                ax_bot.set_ylabel(f'{scalar_labels.get(scalar, scalar)}\n(Group Mean ± 95% CI)', fontsize=10)
                ax_bot.legend(fontsize=8, loc='upper left')
            ax_bot.set_xlabel('Node')

        plt.tight_layout()
        out_path = f'{out_dir}/{tract}_{scalar}_cohort_comparison.png'
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")

# Plot 2: Site-level comparison for key tracts
print("\nGenerating site comparison plots...")
key_tracts = ['C_L', 'C_R', 'F_L', 'F_R', 'UF_L', 'UF_R']
key_tracts = [t for t in key_tracts if t in tracts_to_plot]

for scalar in scalars_to_plot:
    for tract in key_tracts:
        fig, axes = plt.subplots(1, len(approach_names), figsize=(5*len(approach_names), 4),
                                 sharey=True)
        fig.suptitle(f'{tract} — {scalar_labels.get(scalar, scalar)} — By Site', fontsize=14, fontweight='bold')

        for col, approach in enumerate(approach_names):
            df = datasets[approach]
            tract_df = df[df['tractID'] == tract].copy()
            tract_df['nodeID'] = pd.to_numeric(tract_df['nodeID'], errors='coerce')
            tract_df[scalar] = pd.to_numeric(tract_df[scalar], errors='coerce')

            ax = axes[col]
            sites = sorted(tract_df['site'].unique())
            cmap = plt.cm.Set1
            site_colors = {s: cmap(i / max(len(sites)-1, 1)) for i, s in enumerate(sites)}

            for site in sites:
                site_data = tract_df[tract_df['site'] == site]
                grouped = site_data.groupby('nodeID')[scalar]
                mean = grouped.mean()
                sem = grouped.sem()
                nodes = mean.index

                ax.plot(nodes, mean.values, color=site_colors[site], linewidth=2, label=site)
                ax.fill_between(nodes, (mean - 1.96*sem).values, (mean + 1.96*sem).values,
                                color=site_colors[site], alpha=0.15)

            ax.set_title(approach, fontsize=12)
            ax.set_xlabel('Node')
            if col == 0:
                ax.set_ylabel(f'{scalar_labels.get(scalar, scalar)}', fontsize=10)
                ax.legend(fontsize=7, loc='upper left')

        plt.tight_layout()
        out_path = f'{out_dir}/{tract}_{scalar}_site_comparison.png'
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")

print(f"\nAll plots saved to {out_dir}/")
print("Done!")
