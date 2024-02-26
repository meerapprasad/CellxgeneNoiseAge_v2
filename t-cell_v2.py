import anndata
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from tabula_senis.plots import plot_heatmap, plot_count_distributions, plot_linear_regression
from tabula_senis.filter import filter_data, filter_ctype, filter_tissue, filter_mean_cols, filter_coeff_var
from tabula_senis.regression import linear_regression_least_squares


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
    pd.DataFrame(adata_tissue.obs['age']).value_counts().loc[sorted_ages_mm].to_csv(os.path.join(output_dir, 'age_counts.csv'))

    # plot the distribution of counts for each age group
    slope, intercept, mean_filter_df = linear_regression_least_squares(adata_tissue, sorted_ages_mm, output_dir)
    top_genes, bot_genes = plot_linear_regression(slope, intercept, mean_filter_df, output_dir, n_plots=100,
                                                  ncols=8)

    for age in [2,3,4,5]:
        try:
            plot_count_distributions(adata_tissue, age, age_dict, top_genes, sorted_ages_mm, output_dir,
                                     fname='increase')
            plot_count_distributions(adata_tissue, age, age_dict, bot_genes, sorted_ages_mm, output_dir,
                                     fname='decrease')
        except Exception as e:
            continue


if __name__ == "__main__":
    # Example parameters
    cell_type = "T cell"
    tissue = "lung"
    # sorted_ages_mm = ['1m', '3m', '18m', '21m', '24m', '30m']
    adata = anndata.read_h5ad('/home/mprasad/Downloads/2d6aefde-d7ab-4cd0-88d7-f4fba6f94504.h5ad')
    age_dict = {0: '1m', 1: '3m', 2: '18m', 3: '21m', 4: '24m', 5: '30m'}

    main(cell_type, tissue, age_dict, adata)
