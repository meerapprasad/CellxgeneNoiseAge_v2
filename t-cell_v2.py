import anndata
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gc
import numpy as np
from tabula_senis.plots import plot_heatmap, plot_count_distributions, plot_linear_regression
from tabula_senis.filter import filter_data, filter_ctype, filter_tissue, filter_age, combine_labels_adata, filter_coeff_var
from tabula_senis.regression import linear_regression_least_squares, linear_regression_on_means

from DAVIDpy.davidpy import DAVID_start, get_clustering

def tissues_to_keep(adata, cell_type):
    """filter tissues that do not have data at 3m, and filter the age_dict for tissues that do not have data for old ages"""
    adata_sub = filter_ctype(adata, cell_type)
    adata_sub = filter_age(adata_sub)
    contingency_table = pd.crosstab(adata_sub.obs['tissue'], adata_sub.obs.age)
    condition1 = contingency_table['3m'] != 0  # 3m is 0
    condition2 = (contingency_table['18m'] != 0) | (contingency_table['21m'] != 0)  # Both 18m and 21m are 0
    condition3 = (contingency_table['24m'] != 0) | (contingency_table['30m'] != 0)  # Both 24m and 30m are 0
    to_keep = condition1 & (condition2 | condition3)
    # todo: which tissues do we need to filter the sorted_ages_mm
    filter_ages = condition2 & condition3
    return to_keep[to_keep].index.values.tolist(), filter_ages[~filter_ages].index.values.tolist()


def main(cell_type, tissue, age_dict, adata, age_map_dict, age_col):
    # Example of using the plot_heatmap function
    output_dir = f'output/{cell_type}_plots/log-reg-on-means'
    os.makedirs(output_dir, exist_ok=True)
    sorted_ages_mm = list(age_dict.values())

    # Filter data
    adata = filter_data(adata)
    # filter by cell type and plot
    adata_sub = filter_ctype(adata, cell_type)
    # dont consider 4w bc developing
    adata_sub = filter_age(adata_sub)
    adata_sub = combine_labels_adata(adata_sub, age_map_dict, age_col)
    # make heatmap for tissue distributions per age per tissue
    # todo: use this to filter the tissues
    heatmap_filepath = os.path.join(output_dir, 'tissue_type_dist_per_age_heatmap_group.png')
    if not os.path.exists(heatmap_filepath):
        contingency_table = pd.crosstab(adata_sub.obs['tissue'], adata_sub.obs[age_col])
        plot_heatmap(contingency_table, sorted_ages_mm, 'tissue type dist per age', heatmap_filepath)

    output_dir += f'/{tissue}'
    os.makedirs(output_dir, exist_ok=True)
    adata_tissue = filter_tissue(adata_sub, tissue)
    pd.DataFrame(adata_tissue.obs[age_col]).value_counts().loc[sorted_ages_mm].to_csv(os.path.join(output_dir, 'age_counts.csv'))
    thresh = 0.005
    # plot the distribution of counts for each age group
    # slope, intercept, mean_filter_df = linear_regression_least_squares(adata_tissue, sorted_ages_mm, age_col, output_dir)
    # todo: dont run these if there are already plots there
    slope, intercept, mean_filter_df = linear_regression_on_means(adata_tissue, sorted_ages_mm, age_col, output_dir)
    top_genes, bot_genes = plot_linear_regression(slope, intercept, mean_filter_df, output_dir, n_plots=100,
                                                  ncols=8, thresh=thresh)

    for age in [1,2]:
        try:
            plot_count_distributions(adata_tissue, age, age_dict, top_genes, sorted_ages_mm, age_col, output_dir,
                                     fname='increase')
            plot_count_distributions(adata_tissue, age, age_dict, bot_genes, sorted_ages_mm, age_col, output_dir,
                                     fname='decrease')
        except Exception as e:
            continue

    del adata
    gc.collect()
    # todo: run david on the top increase/decrease genes
    increase_genes = pd.read_csv(f'{output_dir}/increase_slope_thresh_{thresh}.csv').sort_values(
        'slope', ascending=False)[:2999].iloc[:, 0]  # .values.tolist()
    decrease_genes = pd.read_csv(f'{output_dir}/decrease_slope_thresh_{thresh}.csv').sort_values(
        'slope', ascending=True)[:2999].iloc[:, 0]  # .values.tolist()
    # todo: maybe run on a different core?
    time.sleep(10)
    output_dir += '/david_results'
    os.makedirs(output_dir, exist_ok=True)
    # todo: save the increase/decrease lists in the david folder, try to input a list with -i
    os.makedirs(os.path.join(output_dir, 'increase'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'decrease'), exist_ok=True)
    increase_genes_list = adata_sub.var[adata_sub.var.feature_name.isin(increase_genes)].index.values.tolist()
    client = DAVID_start(increase_genes_list, save_path=os.path.join(output_dir, 'increase'), email='mprasad@caltech.edu')
    cluster_df = get_clustering(client)
    cluster_df.to_csv(os.path.join(output_dir, 'increase', 'david_results.csv'), index=False)
    # david requires 10 sec interval between hits: https://david.ncifcrf.gov/content.jsp?file=DAVID_API.html
    time.sleep(10)
    decrease_genes_list = adata_sub.var[adata_sub.var.feature_name.isin(decrease_genes)].index.values.tolist()
    client = DAVID_start(decrease_genes_list, save_path=os.path.join(output_dir, 'decrease'), email='mprasad@caltech.edu')
    cluster_df = get_clustering(client)
    cluster_df.to_csv(os.path.join(output_dir, 'decrease', 'david_results.csv'), index=False)


# todo: lots of bugs, because some of the ages are not in the adata
if __name__ == "__main__":
    # Example parameters
    cell_type = "T cell"
    # tissue = "bone marrow"
    # sorted_ages_mm = ['1m', '3m', '18m', '21m', '24m', '30m']
    adata = anndata.read_h5ad('/home/mprasad/Downloads/2d6aefde-d7ab-4cd0-88d7-f4fba6f94504.h5ad')
    # age_dict = {0: '1m', 1: '3m', 2: '18m', 3: '21m', 4: '24m', 5: '30m'}
    # age_dict = {0: '3 month-old stage', 1:'18 month-old stage', 2:'20 month-old stage and over'}
    # age_col = 'development_stage'
    tissue_list, reformat_age_tissue = tissues_to_keep(adata, cell_type)
    for tissue in tissue_list[-2:]:
        if tissue in reformat_age_tissue:
            print(tissue)
            age_map_dict = {'3m': '3m', '18m': '18-21m', '21m': '18-21m'}
            age_col = 'combined_ages'
            age_dict = {0: '3m', 1: '18-21m'}
            main(cell_type, tissue, age_dict, adata, age_map_dict, age_col)
        else:
            age_map_dict = {'3m': '3m', '18m': '18-21m', '21m': '18-21m', '24m': '24-30m', '30m': '24-30m'}
            age_col = 'combined_ages'
            age_dict = {0: '3m', 1: '18-21m', 2: '24-30m'}
            main(cell_type, tissue, age_dict, adata, age_map_dict, age_col)


# todo: ignore tissues that dont have data for 3m, and filter the age_dict for tissues that dont have data for old ages
# todo: instead of saving top and bottom genes only, also save genes where slope is not changing

