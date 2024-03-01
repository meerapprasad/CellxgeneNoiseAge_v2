import os
import pandas as pd
import numpy as np
import sys
import inspect
import re
import configparser
from suds.client import Client

# https://github.com/narek01/DAVIDpy/blob/main/davidpy/davidpy.py
# https://david.ncifcrf.gov/content.jsp?file=WS.html
def check_converter():
    # Checks for existing of file for converting Gene names to Ensembl IDs and reverse

    dir_path = os.path.split(os.path.abspath(inspect.getsourcefile(check_converter)))[0]
    converter_filename = ".ensembl_to_hgnc_human.txt"
    if os.path.isfile(dir_path + '/' + converter_filename):
        print("Ensembl file found")
        return pd.read_csv(dir_path + '/' + converter_filename, sep='\t')
    else:
        try:
            import biomart, sys
            server = biomart.BiomartServer("http://uswest.ensembl.org/biomart")
            # ensembl = server.datasets['hsapiens_gene_ensembl']
            ensembl = server.datasets['mmusculus_gene_ensembl']
            r = ensembl.search({'attributes': ['external_gene_name', 'ensembl_gene_id']})
            convert_df = pd.read_csv(r.url, sep='\t', names=["Gene name", "Gene stable ID"]).dropna()
            convert_df.to_csv(dir_path + '/' + converter_filename, sep='\t', index=False)

            if sys.platform.startswith("win"):  # Hiding file in Windows
                try:
                    import subprocess
                    subprocess.call(['attrib', '+h', dir_path + '/' + converter_filename])
                except:
                    pass

            print("Ensembl file downloaded from Biomart")
            return convert_df
        except:
            print("Couldn't find or download Ensembl file")
            return False


def converter(genes_list_or_path, save_path=None, reverse=False):
    import re

    # Checking if input is a file and trying to open it
    if isinstance(genes_list_or_path, str):
        # todo: the if may not work
        if os.path.isfile(genes_list_or_path):
            with open(genes_list_or_path) as input_file:
                genes_string = input_file.read()
                # Splitting string to list of genes
                genes_list = re.findall(r"[\w\-\_']+", genes_string)
    elif isinstance(genes_list_or_path, list):
        genes_list = genes_list_or_path



    global ensembl_table
    global temporary_conv_df

    # Checking if dataframe already exsists
    try:
        ensembl_table
    except:
        ensembl_table = check_converter()

    # Defining the source of table with gene IDs
    # If we cannot gain a table with ALL genes
    # (e.g. no rights to write file)...
    if type(ensembl_table) != bool:
        temporary_conv_df = ensembl_table
    # ... trying to create the table with OUR genes
    else:
        # if the conversion file already exists, load it instead of re-running the conversion
        try:
            temporary_conv_df = pd.read_csv(f'{save_path}/ensembl_id_conversion.csv')
        except:
            print("Connected to Biomart")
            import biomart
            server = biomart.BiomartServer("http://uswest.ensembl.org/biomart")
            # ensembl = server.datasets['hsapiens_gene_ensembl']
            ensembl = server.datasets['mmusculus_gene_ensembl']
            # iteratively get the genes to avoid query size error
            temporary_conv_df_lst = []
            for i in np.arange(0, len(genes_list), 50):
                r = ensembl.search({'filters': {'external_gene_name': genes_list[i:i + 50]},
                                    'attributes': ['external_gene_name', 'ensembl_gene_id']})
                temporary_conv_df_lst.append(
                    pd.read_csv(r.url, sep='\t', names=["Gene name", "Gene stable ID"]).dropna())
            temporary_conv_df = pd.concat(temporary_conv_df_lst)
            if save_path is not None:
                temporary_conv_df.to_csv(f'{save_path}/ensembl_id_conversion.csv', index=False)

    if not reverse:
        genes_list_conv = temporary_conv_df['Gene stable ID'].loc[
            temporary_conv_df['Gene name'].isin(genes_list)].to_list()
    else:
        genes_list_conv = temporary_conv_df['Gene name'].loc[
            temporary_conv_df['Gene stable ID'].isin(genes_list)].to_list()
    # todo: save these to access later
    return list(set(genes_list_conv))


def set_config(**kwargs):
    import configparser

    if not sys.platform.startswith('linux'):
        print("Your OS is not Linux. If you want to change defaults, edit module's .py file.")
        return
    else:
        path = os.path.expanduser("~") + '/.config/DAVID.ini'
        config = configparser.ConfigParser()
        if os.path.exists(path):
            config.read_file(open(path))
        else:
            config["DEFAULT"] = {"email": "",
                                 "threshold": 0.1,
                                 "count": 2,
                                 "overlap": 3,
                                 "initialSeed": 3,
                                 "finalSeed": 3,
                                 "linkage": 0.5,
                                 "kappa": 50}
            if not "email" in kwargs.keys():
                email = str(input("Please, insert your email: "))
                config["DEFAULT"]["email"] = email

    # Changing parameters
    if kwargs:
        for key in config["DEFAULT"]:
            if key in kwargs.keys():
                config["DEFAULT"][key] = str(kwargs[key])
                del kwargs[key]
                # print(key+":", config["DEFAULT"][key])
    with open(path, 'w') as configfile:
        config.write(configfile)
    if kwargs.keys():
        print("Wrong parameter(s):", ", ".join(kwargs.keys()))


def check_config():
    import configparser
    path = ''
    if sys.platform.startswith('linux'):
        path = os.path.expanduser("~") + '/.config/DAVID.ini'
        if not os.path.exists(path):
            set_config()
        config = configparser.ConfigParser()
        config.read_file(open(path))
        config = dict(config["DEFAULT"])
    else:
        # Change this part if you use OS except of Linux
        config = {"email": "",
                  "threshold": 0.1,
                  "count": 2,
                  "overlap": 3,
                  "initialSeed": 3,
                  "finalSeed": 3,
                  "linkage": 0.5,
                  "kappa": 50}
        if not config['email']:
            email = str(
                input("Please, insert your email and register at http://david-d.ncifcrf.gov/webservice/register.htm\n"))
            config["email"] = email
    return config

# used to take input_list_path
def DAVID_start(input_list, save_path=None, input_bg_path=None, email=None):
    # Reading configs
    global DAVID_conf
    DAVID_conf = check_config()
    if not email:
        email = DAVID_conf["email"]

    # Reading file or string with genes and converting to Ensembl format
    # input_list_ids = ",".join(converter(input_list_path, save_path))
    input_list_ids= ",".join(list(set(input_list)))
    # todo: feed in Ensembl IDs
    print('List loaded')
    if input_bg_path:
        input_bg_ids = ",".join(converter(input_bg_path))
        print('Using file background')
    else:
        print('Using default background')

    # Establishing connection and trying to authenticate
    client = Client(
        'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService?'
        'wsdl')
    client.wsdl.services[0].setlocation(
        'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.'
        'DAVIDWebServiceHttpSoap11Endpoint/')
    auth = client.service.authenticate(email)
    if "Failed" in auth:
        auth = "Failed. For user registration, go to http://david-d.ncifcrf.gov/webservice/register.htm"
    print('User Authentication:', auth)

    # Sending list of genes and background to the server
    print('Percentage mapped(list):', client.service.addList(
        input_list_ids, 'ENSEMBL_GENE_ID', 'LIST1', 0))
    if input_bg_path:
        print('Percentage mapped(background):', client.service.addList(
            input_bg_ids, 'ENSEMBL_GENE_ID', 'BACKGROUND1', 1))
    return client


# https://david-d.ncifcrf.gov/content.jsp?file=DAVID_WebService.html
def get_chart(DAVID, threshold=None, count=None):
    """This is exactly the Functional Annotation Chart"""
    if not threshold:
        threshold = DAVID_conf["threshold"]
    if not count:
        count = DAVID_conf["count"]
    chart = DAVID.service.getChartReport(threshold, count)
    # chart = DAVID.service.getGeneClusterReport(DAVID_conf['overlap'], DAVID_conf['initialseed'],
    #                                            DAVID_conf['finalseed'], DAVID_conf['linkage'], DAVID_conf['kappa'])
    df = pd.DataFrame.from_dict(chart)
    df.columns = [i[0] for i in df.loc[0]]
    df = df.applymap(lambda x: x[1])
    new_IDs = []
    for ID in df['geneIds']:
        new_IDs.append(", ".join(converter(ID, reverse=True)))
    df['geneIds'] = new_IDs
    cols = ["categoryName", "termName", "geneIds", "percent", "ease", "benjamini"]
    df = df[cols].sort_values("benjamini")
    df.columns = ["Category", "Term", "Genes", "Percent", "P-Value", "Benjamini"]
    return df


def get_clustering(DAVID, overlap=None, initialSeed=None, finalSeed=None, linkage=None, kappa=None):
    if not overlap:
        overlap = DAVID_conf['overlap']
    if not initialSeed:
        initialSeed = DAVID_conf['initialseed']
    if not finalSeed:
        finalSeed = DAVID_conf['finalseed']
    if not linkage:
        linkage = DAVID_conf['linkage']
    if not kappa:
        kappa = DAVID_conf['kappa']
    # chart = DAVID.service.getGeneClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)
    chart = DAVID.service.getTermClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)
    # todo: make this a dictionary and combine from that ...
    # df = pd.DataFrame.from_dict(chart)
    info_ids = ['termName', 'foldEnrichment', 'geneIds', 'benjamini', 'bonferroni', 'EASEBonferroni', 'fisher', 'rfdr']
    # chart[0]['simpleChartRecords'][0]['foldEnrichment']
    rows_lst = []
    for i in range(len(chart)):
        for x in range(len(chart[i]['simpleChartRecords'])):
            row = []
            for j in info_ids:
                if j == 'geneIds':
                    row.append([chart[i]['simpleChartRecords'][x][j]])
                elif j == 'termName' and ';' in chart[i]['simpleChartRecords'][x][j]:
                     row.append(','.join(chart[i]['simpleChartRecords'][x][j].rsplit(';')))
                else:
                    row.append(chart[i]['simpleChartRecords'][x][j])
            row.append(i)
            rows_lst.append(row)

    df = pd.DataFrame(rows_lst, columns=info_ids+['cluster_id'])
    return df


def main():
    import argparse
    parser = argparse.ArgumentParser(description='DAVID Functional Annotation Chart retrieving script')
    parser.add_argument('-i', type=str, help='Input list of genes (file or "GENE1,GENE2,ETC")', required=True)
    parser.add_argument('--bg', type=str, help='Background (file or "GENE1,GENE2,ETC")')
    parser.add_argument('--tsv', help='Save to tsv', action="store_true")
    parser.add_argument('--csv', help='Save to csv', action="store_true")
    parser.add_argument('--full', help='Show all columns of the table', action="store_true")
    parser.add_argument('--save_path', type=str, help='Path to save the file')

    args = parser.parse_args()

    if not args.i:
        DAVID = DAVID_start(args.i)
    else:
        DAVID = DAVID_start(args.i, args.bg)

    df = get_chart(DAVID)
    cluster_df = get_clustering(DAVID)
    cluster_df.to_csv(f"{args.save_path}/functional_annotation_cluster.csv", index=False)

    print("First 20 terms:")
    if args.full:
        pd.options.display.max_columns = None
        print(df.head(20))
    else:
        print(df.head(20)[['Term', 'Percent', 'Benjamini']])

    if args.tsv:
        df.to_csv("DAVID_chart.tsv", sep='\t', index=False)
    elif args.csv:
        df.to_csv("DAVID_chart.csv", index=False)


if __name__ == '__main__':
    main()