# import requests
# from bs4 import BeautifulSoup
# import pandas as pd
#
# # URL for the page where you submit your gene list
# url = 'https://david.ncifcrf.gov/summary.jsp'
#
# # Example gene list - replace with your actual list
# gene_list = 'BRCA1, BRCA2, TP53'
#
# # Form data to mimic what would be sent by the webpage form
# # Adjust the field names and values based on actual form requirements
# form_data = {
#     'geneList': gene_list,  # Field for pasting the gene list
#     'identifier': 'OFFICIAL_GENE_SYMBOL',  # Hypothetical identifier selection
#     'listType': 'Gene List',  # Assuming you're submitting a gene list, not background
#     'species': 'Mus musculus',
#     'submit': 'Submit'  # You might or might not need to mimic a submit button press
# }
#
# # Make the POST request
# response = requests.post(url, data=form_data)
#
# soup = BeautifulSoup(response.text, 'html.parser')


# (pd.read_csv('output/T cell_plots/bone marrow/decrease_slope_thresh_0.005.csv').iloc[:,0]
#  .to_csv('output/T cell_plots/bone marrow/gene_list.csv', index=False))

import davidpy
client = davidpy.DAVID_start('output/T cell_plots/bone marrow/gene_list.csv', email='mprasad@caltech.edu')
df = davidpy.get_chart(client)
print('done')




## TODO: hw calculations
import numpy as np
# sigma = 5.67e-8
# a_i = (4.2/300)**4
# a_ii = (300**4 + 4.2**4)/(2*77**4)
# a_iii = sigma* (300**4 - 77**4)
#
# b_i = ((77**4 + 4.2**4)/2)**(1/4)
# b_ii = sigma*(77**4-b_i**4)