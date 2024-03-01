from DAVIDpy.davidpy import DAVID_start, get_chart, get_clustering
import pandas as pd

gene_list_file = 'output/T cell_plots/log-reg-least-squares/bone marrow/gene_list.csv'
gene_list = pd.read_csv(gene_list_file).feature_name.values.tolist()
save_path = 'output/T cell_plots/log-reg-least-squares/bone marrow'
client = DAVID_start(gene_list, email='mprasad@caltech.edu')
# df = davidpy.get_chart(client)
cluster_df = get_clustering(client)

print('done')




# hw calculations, lol
import numpy as np
# sigma = 5.67e-8
# a_i = (4.2/300)**4
# a_ii = (300**4 + 4.2**4)/(2*77**4)
# a_iii = sigma* (300**4 - 77**4)
#
# b_i = ((77**4 + 4.2**4)/2)**(1/4)
# b_ii = sigma*(77**4-b_i**4)