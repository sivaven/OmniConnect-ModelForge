################################################################################
# File: figs_aba_scale.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Generates Allen Brain Atlas (ABA) scaling validation figures.
################################################################################

import cross_species_nets as cn
import stat_utils as su
import vis_utils as vu

import numpy as np
import pandas as pd
import networkx as nx
from numpy import inf
import matplotlib.pyplot as plt

_dir='data/'

# Fig 3d
mat_dmri, mat_tracer, gm_labels, mat_count_dmri = cn.get_cmats_aba_scale(root_dir=_dir, atlas='ABAminor')
vols=[]    
vol_white=pd.read_table(_dir+'mouse_wmstats.txt', index_col=0).T['volume (mm^3)'].sum()    
vol_white=vol_white**(1/3)    
mat_dmri=mat_dmri/vol_white    

# Global costs, Null comparisons
spl_dmri, spl_dmri_tracer, spls_rand_asymm_dmri, spls_rand_shuffled_dmri, spls_rand_allshuffle = \
    cn.get_avg_spls(mat_dmri[:144, :144], mat_tracer[:144, :144], n_trials=10, print_progress=True) #set to 10000

# path rankings
mat_asym=1/(mat_tracer)
mat_asym[mat_asym==inf]=0
mat_asym=np.nan_to_num(mat_asym, 0)
g=nx.from_numpy_array(mat_dmri[:144, :144]*mat_asym[:144, :144], create_using=nx.DiGraph) 
spls_dict={}
mat_inv_dmri_count=1/(mat_count_dmri)
mat_inv_dmri_count[mat_inv_dmri_count==inf]=0
mat_inv_dmri_count=np.nan_to_num(mat_inv_dmri_count, 0)                                    
g_sym=nx.from_numpy_array(mat_dmri[:144, :144]*mat_inv_dmri_count, create_using=nx.DiGraph)
spls_dict_sym={}

def main():
    np.fill_diagonal(mat_dmri, 0)
    np.fill_diagonal(mat_count_dmri, 0)
    np.fill_diagonal(mat_tracer, 0)

    for src_id in range(len(gm_labels[:144])):
        try:
            dict_=nx.shortest_path_length(g, source=src_id, weight=cn.weighted_shortest_path)
            for tar_id in range(len(gm_labels[:144])):
                if dict_[tar_id]>0:
                    key=gm_labels[src_id]+'__'+gm_labels[tar_id]
                    spls_dict[key]=dict_[tar_id]
        except nx.NetworkXNoPath:
            continue
    for src_id in range(len(gm_labels[:144])):
        try:
            dict_=nx.shortest_path_length(g_sym, source=src_id, weight=cn.weighted_shortest_path)
            for tar_id in range(len(gm_labels[:144])):
                if dict_[tar_id]>0:
                    key=gm_labels[src_id]+'__'+gm_labels[tar_id]
                    spls_dict_sym[key]=dict_[tar_id]
        except nx.NetworkXNoPath:
            continue

    cost_dict_asym = dict(sorted(spls_dict.items(), key=lambda item: item[1], reverse=True))
    cost_dict_sym  = dict(sorted(spls_dict_sym.items(), key=lambda item: item[1], reverse=True))

    vu.plot_cost_dict_aba_scale([cost_dict_asym, cost_dict_sym])

    # Fig 2d heatmaps
    fig=plt.figure(figsize=(7,3.5))
    plt.subplot(121)
    plt.imshow((mat_count_dmri), vmin=0, vmax=2000)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.gca().set_xticks([25, 50, 75, 100, 125])
    plt.gca().set_yticks([25, 50, 75, 100, 125])
    plt.title('dMRI Tractography\n (Streamline counts)')
    plt.xlabel('ABA ROIs'); plt.ylabel('ABA ROIs')

    plt.subplot(122)
    mat_tr = mat_tracer[:144, :144].copy()
    np.fill_diagonal(mat_tr, 0)
    plt.imshow(mat_tr, vmin=0, vmax=.1)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.gca().set_xticks([25, 50, 75, 100, 125])
    plt.gca().set_yticks([25, 50, 75, 100, 125])
    plt.title('Viral Tracing\n (Projection densities)')
    plt.xlabel('ABA ROIs')
    plt.ylabel('ABA ROIs')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
