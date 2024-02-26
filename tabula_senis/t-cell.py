import anndata
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
# from plotting.colors import godsnot_64
# from matplotlib.colors import to_rgb
# import matplotlib.patches as mpatches
from cluster_methods.hierarchical import hierarchical_cluster
from tabula_senis.plots import (plot_heatmap, plot_count_distributions, plot_linear_regression,
                                plot_count_dist_noisy_gene_per_age)
from tabula_senis.filter import filter_data, filter_ctype, filter_tissue, filter_mean_cols, filter_coeff_var

def compute_log_mean_and_var(adata_tissue, sorted_ages_mm):
    # fano factor: var/ mean
    log_mean_df = pd.DataFrame(np.zeros((len(adata_tissue.obs.age.unique().tolist()), len(adata_tissue.var_names))),
                               index=sorted_ages_mm, columns=adata_tissue.var['feature_name'])
    for val in sorted_ages_mm:
        log_mean_df.loc[val, :] = np.log1p(adata_tissue[adata_tissue.obs.age == val].X.todense()).mean(axis=0)

    log_var_df = pd.DataFrame(np.zeros((len(adata_tissue.obs.age.unique().tolist()), len(adata_tissue.var_names))),
                              index=sorted_ages_mm, columns=adata_tissue.var['feature_name'])
    for val in sorted_ages_mm:
        log_data = np.log1p(adata_tissue[adata_tissue.obs.age == val].X.todense())
        # Compute squared deviations from the mean
        squared_deviations = (np.array(log_data) - np.array(log_mean_df.loc[val, :])) ** 2
        # Sum squared deviations and divide by the number of cells to get variance
        variance = np.sum(squared_deviations, axis=0) / (log_data.shape[0] - 1)
        # Assign variance to the corresponding row in log_var_df
        log_var_df.loc[val, :] = variance

    return log_mean_df, log_var_df

def compute_and_save_statistical_measures(adata, sorted_ages, output_dir):
    log_mean_df, log_var_df = compute_log_mean_and_var(adata, sorted_ages)
    ff_df = (log_var_df / log_mean_df).fillna(0)
    cv_df = (np.sqrt(log_var_df) / log_mean_df).fillna(0)
    ff_df.to_csv(os.path.join(output_dir, 'fano_factor_per_gene_per_age.csv'), index=False)
    cv_df.to_csv(os.path.join(output_dir, 'coeffofvar_per_gene_per_age.csv'), index=False)
    return ff_df, cv_df


def hierarchical_clustering_and_plotting(df, sorted_ages, output_dir, title):
    # Assume hierarchical_cluster is a predefined function
    clusters = hierarchical_cluster(df.T.values)
    np.save(os.path.join(output_dir, 'cv_hclust.npy'), np.stack([df.columns.values, clusters]).T)
    plt.imshow(df.values[:, clusters], aspect='auto', cmap='coolwarm')
    plt.yticks(ticks=np.arange(df.shape[0]), labels=sorted_ages)
    plt.colorbar()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cv_per_gene_per_age.png'), dpi=300)
    plt.clf()
    return clusters

def main(cell_type, tissue, age_dict, adata):
    # Example of using the plot_heatmap function
    output_dir = f'output/{cell_type}_plots'
    os.makedirs(output_dir, exist_ok=True)
    sorted_ages_mm = list(age_dict.values())

    # Filter data
    adata = filter_data(adata)
    # todo: filter by cell type and plot
    adata_sub = filter_ctype(adata, cell_type)
    # make heatmap for tissue distributions per age per tissue
    heatmap_filepath = os.path.join(output_dir, 'tissue_type_dist_per_age_heatmap.png')
    if not os.path.exists(heatmap_filepath):
        contingency_table = pd.crosstab(adata_sub.obs['tissue'], adata_sub.obs['age'])
        plot_heatmap(contingency_table, sorted_ages_mm, 'tissue type dist per age', heatmap_filepath)

    output_dir += f'/{tissue}'
    os.makedirs(output_dir, exist_ok=True)
    adata_tissue = filter_tissue(adata_sub, tissue)
    ff_df, cv_df = compute_and_save_statistical_measures(adata_tissue, sorted_ages_mm, output_dir)
    df_filtered = filter_coeff_var(cv_df)

    if not os.path.exists(os.path.join(output_dir, 'cv_hclust.npy')):
        clusters = hierarchical_clustering_and_plotting(df_filtered, sorted_ages_mm, output_dir,
                                                        'coeff of var gene per age')
    else:
        clusters = np.load(os.path.join(output_dir, 'cv_hclust.npy'), allow_pickle=True)[:, 1]
    # todo: plot the distribution of counts for each age group

    for age in [2,3,4,5]:
        try:
            plot_count_dist_noisy_gene_per_age(df_filtered, adata_tissue, sorted_ages_mm, clusters, age, age_dict,
                                               output_dir)
        except Exception as e:
            continue


if __name__ == "__main__":
    # Example parameters
    cell_type = "T cell"
    tissue = "kidney"
    # sorted_ages_mm = ['1m', '3m', '18m', '21m', '24m', '30m']
    adata = anndata.read_h5ad('/home/mprasad/Downloads/2d6aefde-d7ab-4cd0-88d7-f4fba6f94504.h5ad')
    age_dict = {0: '1m', 1: '3m', 2: '18m', 3: '21m', 4: '24m', 5: '30m'}

    main(cell_type, tissue, age_dict, adata)
