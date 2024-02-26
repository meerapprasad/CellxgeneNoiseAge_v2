import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import pandas as pd

def plot_heatmap(data, columns, title, filepath, figsize=(8, 6)):
    plt.figure(figsize=figsize)
    heatmap_data = np.log1p(data.reindex(columns=columns))
    sns.heatmap(heatmap_data, annot=False, cbar=True, fmt=".2f")
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.clf()


def plot_count_distributions(adata_tissue, age, age_dict, gene_names, sorted_ages_mm, output_dir, fname='increase'):
    n_genes = gene_names.shape[0]  # Number of genes to plot
    n_cols = 6  # Number of columns in the subplot # grid
    n_rows = (n_genes + n_cols - 1) // n_cols  # Calculate rows needed
    os.makedirs(os.path.join(output_dir, 'count_dist_noisy_genes'), exist_ok=True)
    plt.figure(figsize=(5 * n_cols, 4 * n_rows))
    # age_groups = adata_tissue.obs['age'].unique()
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'yellow']  # [:len(age_groups)]
    # make a folder to save each one
    os.makedirs(os.path.join(output_dir, 'count_dist_noisy_genes'), exist_ok=True)
    sub_ages = ['1m', '3m', sorted_ages_mm[age]]
    sub_colors = [colors[0], colors[1], colors[age]]
    # plot only 1m, 3m and final age group
    # todo: normalize based on gene counts ?
    for i, gene_idx in enumerate(np.where(adata_tissue.var['feature_name'].values.isin(gene_names))[0]):
        plt.subplot(n_rows, n_cols, i + 1)
        # for age_group, color in zip(sorted_ages_mm, colors):
        for age_group, color in zip(sub_ages, sub_colors):
            # Select the cells for a given age group and gene
            vals = np.log1p(
                np.array(adata_tissue[adata_tissue.obs['age'] == age_group, gene_idx].X.todense()).flatten())
            # plt.hist(vals, bins=15, alpha=0.7, label=str(age_group), color=color, log=True, density=True)
            weights = np.ones_like(vals) / len(vals)
            # Plot the histogram with weights
            plt.hist(vals, bins=15, alpha=0.7, label=str(age_group), color=color, weights=weights)

        plt.title(f'{adata_tissue.var["feature_name"].values[gene_idx]}')
        plt.xlabel('log(counts + 1)')
        plt.ylabel('frequency')
        if i == 0:  # Add legend only to the first subplot for clarity
            plt.legend(title='Age Group')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/count_dist_noisy_genes/age_{age_dict[age]}_sub_norm_{fname}.png', dpi=300)
    plt.clf()


# todo: plot all genes with non-zero slope
def plot_linear_regression(slope, intercept, mean_filter_df, output_dir, n_plots=100, ncols=8, thresh=.005):
    # Sort the slopes and intercepts
    slope_sorted_idx = np.argsort(slope)
    intercept_sorted_idx = np.argsort(intercept)
    # Get the top n_plots genes with the highest and lowest slopes and intercepts
    top_n_slope = slope_sorted_idx[-n_plots:]
    bottom_n_slope = slope_sorted_idx[:n_plots]

    # todo: save the top and bot genes above and below a threshold
    pd.DataFrame(slope[(slope > thresh)]).to_csv(f'{output_dir}/increase_slope_thresh_{thresh}.csv')
    pd.DataFrame(slope[(slope < -thresh)]).to_csv(f'{output_dir}/decrease_slope_thresh_{thresh}.csv')
    # Plot the genes with the top n_plots slopes
    n_rows = int(np.ceil(len(top_n_slope) / ncols))
    time_points = np.array(mean_filter_df.index.str.extract('(\d+)m')[0], dtype=int)

    plt.figure(figsize=(5 * ncols, 4 * n_rows))  # Adjust figsize as needed
    for i, idx in enumerate(top_n_slope):
        plt.subplot(n_rows, ncols, i + 1)  # Update subplot to use n_rows and 8 columns
        plt.plot(time_points, mean_filter_df.iloc[:, idx], 'o')  # Original data points
        plt.plot(np.arange(time_points.min(), time_points.max()),
                 slope[idx] * np.arange(time_points.min(), time_points.max()) + intercept[idx],
                 '-')  # Regression line
        plt.title(f'{mean_filter_df.columns[idx]}', fontsize=10)  # You might need to adjust font size for clarity
        plt.xlabel('Age (months)', fontsize=8)
        plt.ylabel('Mean log(counts + 1)', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'increase_slope_genes_{n_plots}.png'), dpi=300)
    plt.clf()
    # Plot the genes with the bottom n_plots slopes
    n_rows = int(np.ceil(len(bottom_n_slope) / ncols))

    plt.figure(figsize=(5 * ncols, 4 * n_rows))  # Adjust figsize as needed
    for i, idx in enumerate(bottom_n_slope):
        plt.subplot(n_rows, ncols, i + 1)  # Update subplot to use n_rows and 8 columns
        plt.plot(time_points, mean_filter_df.iloc[:, idx], 'o')  # Original data points
        plt.plot(np.arange(time_points.min(), time_points.max()),
                 slope[idx] * np.arange(time_points.min(), time_points.max()) + intercept[idx],
                 '-')  # Regression line
        plt.title(f'{mean_filter_df.columns[idx]}', fontsize=10)  # Adjust title font size
        plt.xlabel('Age (months)', fontsize=8)
        plt.ylabel('Mean log(counts + 1)', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'decrease_slope_genes_{n_plots}.png'), dpi=300)
    plt.clf()  # Clear the figure
    return (mean_filter_df.iloc[:, top_n_slope].columns.to_numpy(),
            mean_filter_df.iloc[:, bottom_n_slope].columns.to_numpy())


def plot_count_dist_noisy_gene_per_age(df_filtered, adata_tissue, sorted_ages_mm, clusters, age, age_dict, output_dir):
    ordered_df = df_filtered.iloc[:, clusters.tolist()]
    bool_idx = adata_tissue.var['feature_name'].values.isin(ordered_df.iloc[:, -500:].columns.to_numpy()[
                                                                np.argmax(ordered_df.iloc[:, -500:], axis=0) == age])
    n_genes = np.sum(bool_idx)  # Number of genes to plot
    n_cols = 6  # Number of columns in the subplot # grid
    n_rows = (n_genes + n_cols - 1) // n_cols  # Calculate rows needed

    plt.figure(figsize=(5 * n_cols, 4 * n_rows))
    # age_groups = adata_tissue.obs['age'].unique()
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'yellow']  # [:len(age_groups)]
    # make a folder to save each one
    os.makedirs(os.path.join(output_dir, 'count_dist_noisy_genes'), exist_ok=True)
    sub_ages = ['1m', '3m', sorted_ages_mm[age]]
    sub_colors = [colors[0], colors[1], colors[age]]
    # plot only 1m, 3m and final age group
    # todo: normalize based on gene counts ?
    for i, gene_idx in enumerate(np.where(bool_idx)[0]):
        plt.subplot(n_rows, n_cols, i + 1)
        # for age_group, color in zip(sorted_ages_mm, colors):
        for age_group, color in zip(sub_ages, sub_colors):
            # Select the cells for a given age group and gene
            vals = np.log1p(
                np.array(adata_tissue[adata_tissue.obs['age'] == age_group, gene_idx].X.todense()).flatten())
            # plt.hist(vals, bins=15, alpha=0.7, label=str(age_group), color=color, log=True, density=True)
            weights = np.ones_like(vals) / len(vals)
            # Plot the histogram with weights
            plt.hist(vals, bins=15, alpha=0.7, label=str(age_group), color=color, weights=weights)

        plt.title(f'{adata_tissue.var["feature_name"].values[gene_idx]}')
        plt.xlabel('log(counts + 1)')
        plt.ylabel('frequency')
        if i == 0:  # Add legend only to the first subplot for clarity
            plt.legend(title='Age Group')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/count_dist_noisy_genes/age_{age_dict[age]}_sub_norm.png', dpi=300)
