#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy
import os, glob
import random

###--------------------------------------------4 layers

adata =  sc.read("/nemo/lab/gandhis/home/users/grantpm/DLPFC_data/DLPFC_Visium/combined_dataframes/nonIF_c2l_anndata_combined_30samples.h5ad")

# Prepare data
subset = adata[adata.obs['bayesSpace_harmony_7'].isin([1,2,3,4,5,6])]   # Excluding WM and meningeal layers

subset.obs['bayesSpace_harmony_7'] = subset.obs['bayesSpace_harmony_7'].replace({   # Renaming to match order of MRI (L1, closer to WM, L4 closer to meninges)
                                                     1 : 'L1',
                                                     2 : 'L4', 
                                                     3 : 'L3', 
                                                     4 : 'L2', 
                                                     5 : 'WM',
                                                     6 : 'L1'})

subset.obs['sample_id'] = subset.obs['sample_id'].astype(str)

df_deconvo = subset.obs[np.isnan(subset.obs['Oligo'])==False]

paths = glob.glob('../../processed_data/00-prepare_GO/01-fetch_genes_for_GO_pathways/*.csv')                         # read gene per pathway data

for m in range(2,4):
    df = pd.read_csv(paths[m], index_col=0)
    pathways = df['GO'].unique()
    categories = pd.DataFrame(paths)[0].str.split('/', expand = True)[5].str.split('.csv', expand = True)[0]
    print(len(pathways))
    for n in range(0,len(pathways)):
        genes = list(df[df['GO']==pathways[n]]['SYMBOL'])
        ensembl = list(subset.var[subset.var['gene_name'].isin(genes)].index)
        print(len(genes))
        for it in range(0,1000):
            universe = list(subset.var.index)
            gene_pool = list(set(universe)-set(ensembl))                           
            interest = subset[:, ensembl].to_df().sum(axis = 1)                    #Take the sum of gene expression for genes of interest                     
            r = random.sample(gene_pool, len(ensembl))
            background = subset[:, r].to_df().sum(axis = 1)                                  #Take sum of expression for background genes
            subset.obs['score'] = (interest - background)/subset.obs['sum_umi']
            df1 = subset.obs[['sample_id','bayesSpace_harmony_7','score']]
            if (it == 0):
                all_scores = subset.obs[['bayesSpace_harmony_7','sample_id','score']]
            else: 
                all_scores = pd.concat([all_scores, subset.obs['score']], axis = 1)

        median  = all_scores.groupby(by = ['sample_id', 'bayesSpace_harmony_7']).sum().median(axis=1).reset_index().rename(columns={0 : 'score'})
        sns.boxplot(data=median, x="bayesSpace_harmony_7", y="score")
        sns.swarmplot(data=median, x="bayesSpace_harmony_7", y="score", color = 'black')
        plt.title(categories[m]+' : '+pathways[n])
        plt.savefig('../../plots/01a-4layers_noendo/05-pathway_scoring/scores_4layers_%s_%s.png' % (categories[m],pathways[n]), bbox_inches='tight')
        plt.show()
        plt.clf()
            
        median.to_csv('../../processed_data/01a-4layers_noendo/05-pathway_scoring/scores_4layers_%s_%s.csv' % (categories[m],pathways[n]))
