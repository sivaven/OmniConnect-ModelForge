################################################################################
# File: fig3b-4a.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Generates Figures 3b and 4a. Performs source/target efficiency analyses across species,
#              applies Kruskal–Wallis tests with multiple comparisons correction,
#              and produces evolutionary efficiency plots.
################################################################################

import cross_species_nets as cn
import stat_utils as su
import vis_utils as vu

import numpy as np
import pandas as pd
import networkx as nx
from numpy import inf

from stat_utils import kruskal
from statsmodels.stats.multitest import multipletests

_dir='data/'

list_of_dicts=[]
mat_dmris=[]

alls_s=[]
allt_s=[]
alls_p=[]
allt_p=[]
src_all=[]
tar_all=[]

def main():
    for idx, species in enumerate(['mouse', 'marmoset', 'rhesus', 'human']):
        print(species)
        mat_dmri, mat_tracer, gm_labels=cn.get_cmats(root_dir=_dir, 
                                                      species=species,
                                                      exclude_indices=[15,21,24,25,26])    
        vol_white=pd.read_table(_dir+species+'_wmstats.txt', index_col=0).T['volume (mm^3)'].sum()
        vol_white=vol_white**(1/3)
        mat_dmri=mat_dmri/vol_white
        mat_dmris.append(mat_dmri)

        np.fill_diagonal(mat_dmri, 0)
        np.fill_diagonal(mat_tracer, 0)

        mat_asym=1/(mat_tracer)
        mat_asym[mat_asym==inf]=0
        mat_asym=np.nan_to_num(mat_asym, 0)

        g=nx.from_numpy_array(mat_dmri*mat_asym, create_using=nx.DiGraph) 
        spls_dict={}
        for src_id in range(len(gm_labels)):
            try:
                dict_=nx.shortest_path_length(g, source=src_id, weight=cn.weighted_shortest_path)
                for tar_id in range(len(gm_labels)):
                    if dict_[tar_id]>0:
                        key=gm_labels[src_id]+'__'+gm_labels[tar_id]
                        spls_dict[key]=dict_[tar_id]
            except nx.NetworkXNoPath as e:
                print(e)
        list_of_dicts.append(spls_dict)

    cost_dict_mouse  = dict(sorted(list_of_dicts[0].items(), key=lambda kv: kv[1], reverse=True))
    cost_dict_marm   = dict(sorted(list_of_dicts[1].items(), key=lambda kv: kv[1], reverse=True))
    cost_dict_rhesus = dict(sorted(list_of_dicts[2].items(), key=lambda kv: kv[1], reverse=True))
    cost_dict_human  = dict(sorted(list_of_dicts[3].items(), key=lambda kv: kv[1], reverse=True))

    n_paths=len(gm_labels)*len(gm_labels) - len(gm_labels)
    ascending_ranks = np.arange(1, n_paths + 1)
    df_mouse    = pd.DataFrame(ascending_ranks[::-1], index=list(cost_dict_mouse.keys()),  columns=['mouse'])
    df_marmoset = pd.DataFrame(ascending_ranks[::-1], index=list(cost_dict_marm.keys()),   columns=['marmoset'])
    df_rhesus   = pd.DataFrame(ascending_ranks[::-1], index=list(cost_dict_rhesus.keys()), columns=['rhesus'])
    df_human    = pd.DataFrame(ascending_ranks[::-1], index=list(cost_dict_human.keys()),  columns=['human'])
    df_evolution = df_mouse.join(df_marmoset).join(df_rhesus).join(df_human)
    df_evolution = 1/df_evolution
    df_evolution = (df_evolution - df_evolution.min()) / (df_evolution.max() - df_evolution.min())
    df_evolution = np.log1p(1000 * df_evolution)

    vu.plot_cost_dict([cost_dict_mouse, cost_dict_marm, cost_dict_rhesus, cost_dict_human])

    # roi-wise stats
    for roi in gm_labels:   
        src_avg=[]
        tar_avg=[]        
        for path in df_evolution.index:
            src=path.split('__')[0]
            tar=path.split('__')[1]        
            if src==roi:            
                src_avg.append(df_evolution.loc[path])
            if tar==roi:            
                tar_avg.append(df_evolution.loc[path])
        src_f=kruskal(pd.DataFrame(src_avg)['mouse'],
                      pd.DataFrame(src_avg)['marmoset'],
                      pd.DataFrame(src_avg)['rhesus'],
                      pd.DataFrame(src_avg)['human'])    
        tar_f=kruskal(pd.DataFrame(tar_avg)['mouse'],
                      pd.DataFrame(tar_avg)['marmoset'],
                      pd.DataFrame(tar_avg)['rhesus'],
                      pd.DataFrame(tar_avg)['human'])   
        alls_s.append(src_f[0]); alls_p.append(src_f[1]); src_all.append(src_avg)
        allt_s.append(tar_f[0]); allt_p.append(tar_f[1]); tar_all.append(tar_avg)

    adjusted_p_values = multipletests(alls_p+allt_p, method='fdr_bh')[1]
    all_adj_s_p = list(adjusted_p_values[:len(gm_labels)])
    all_adj_t_p = list(adjusted_p_values[len(gm_labels):])

    vu.plot_ev(src_all, all_adj_s_p, alls_s, 'Source', gm_labels=gm_labels,
               species_names=['Mouse','Marmoset','Rhesus','Human'])
    vu.plot_ev(tar_all, all_adj_t_p, allt_s, 'Target', gm_labels=gm_labels,
               species_names=['Mouse','Marmoset','Rhesus','Human'])

if __name__ == '__main__':
    main()
