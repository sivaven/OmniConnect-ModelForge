################################################################################
# File: cross_species_nets.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Core library for constructing directed connectomes across species.
#              Implements path efficiency metric, weighted shortest-path algorithms,
#              and network utilities for comparative analyses.
################################################################################
import networkx as nx
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import collections
import pandas as pd
from numpy import inf

def get_cmats(root_dir, 
              species, 
              exclude_indices=[]):    
    """Auto-added docstring. See manuscript methods for details."""
    gm_atlas=root_dir+'moused_'+species+'_CHA.nii.gz'
    gm_labels=pd.read_table(gm_atlas.replace('nii.gz', 'txt'), 
                            delim_whitespace=True, 
                            header=None,
                           index_col=0)[1].values    
    mat_dmri=scipy.io.loadmat(root_dir+'moused_'+species+'_CHA_mean_length_dmri.mat')['connectivity']   
    mat_tracer=pd.read_csv(root_dir+'mouse_CHA__tracer.csv', index_col=0).values
    
    mat_tracer=np.delete(mat_tracer, exclude_indices, axis=0)
    mat_tracer=np.delete(mat_tracer, exclude_indices, axis=1)        
    gm_labels=np.delete(gm_labels, exclude_indices)
    mat_dmri=np.delete(mat_dmri, exclude_indices, axis=0)
    mat_dmri=np.delete(mat_dmri, exclude_indices, axis=1)    
    
    return mat_dmri, mat_tracer, gm_labels

def get_cmats_aba_scale(root_dir, atlas):    
    """Auto-added docstring. See manuscript methods for details."""
    gm_atlas=root_dir+atlas+'.nii.gz'
    gm_labels=pd.read_table(gm_atlas.replace('nii.gz', 'txt'), 
                            delim_whitespace=True, 
                            header=None,
                           index_col=0)[1].values    
    mat_dmri=scipy.io.loadmat(root_dir+atlas+'_mean_length_dmri.mat')['connectivity']   
    mat_count_dmri=scipy.io.loadmat(root_dir+atlas+'_count_dmri.mat')['connectivity']   
    mat_tracer=pd.read_csv(root_dir+atlas+'__tracer.csv', index_col=0).values
      
    return mat_dmri, mat_tracer, gm_labels, mat_count_dmri
    
    
def pathlist2edges(source, lst, gmlabels):
    """Auto-added docstring. See manuscript methods for details."""
    edges=[]
    if len(lst)>0:
        edge=gmlabels[lst[0]]
        for idx, item in enumerate(lst[1:]):
            edges.append(edge+'___'+gmlabels[item])
            edge=gmlabels[item]
        
    return edges
    
def weighted_shortest_path(u, v, edge_attributes):      
    """Auto-added docstring. See manuscript methods for details."""
    return edge_attributes['weight']
#
#
def directed_mspl(g, nodeid=-1):    
    """Auto-added docstring. See manuscript methods for details."""
    if nodeid==-1:
        n_nodes=g.number_of_nodes()
        sum_of_spls=0
        splss=[]
        for spl_dict in nx.shortest_path_length(g, weight=weighted_shortest_path):
            spls=[spl for spl in spl_dict[1].values() if spl>0]
            splss.append(spls)
        splss = [
        ii
        for i in splss
        for ii in i
        ]
        #plt.boxplot(splss)
        #plt.show()
        return splss 
    else:
        lst=[np.inf]
        try:
            dict_=nx.shortest_path_length(g, source=nodeid, weight=weighted_shortest_path)
            lst=list(dict_.values())
        except nx.NetworkXNoPath as e:
            print(e)
        return np.sum(lst)
    
def shuffle_symmetric_matrix(matrix):
    """Auto-added docstring. See manuscript methods for details."""
    assert np.allclose(matrix, matrix.T), "Matrix is not symmetric"
    perm = np.random.permutation(matrix.shape[0])    
    shuffled_matrix = matrix[perm][:, perm]    
    return shuffled_matrix

def symmetrize_tracer_mat(tracer_mat):
    """Auto-added docstring. See manuscript methods for details."""
    tracer_mat_copy=tracer_mat.copy()
    for i in range(len(tracer_mat_copy)):
        for j in range(len(tracer_mat_copy[i])):
            if j<=i:
                continue
            tracer_mat_copy[i][j]=(tracer_mat_copy[i][j]+tracer_mat_copy[j][i])/2
            tracer_mat_copy[j][i]=tracer_mat_copy[i][j]
    return tracer_mat_copy

def shuffle_tracer_preserve_row_idx(mat, seed=None):
    """Auto-added docstring. See manuscript methods for details."""
    rng = np.random.default_rng(seed)
    out = mat.copy()
    for i in range(out.shape[0]):
        nz = np.flatnonzero(out[i, :])
        if nz.size > 1:
            perm = rng.permutation(nz.size)
            out[i, nz] = out[i, nz][perm]
    return out
    
def get_avg_spls(mat_dmri, 
                 mat_tracer, 
                 nodeid=-1, 
                 n_trials=100, 
                 print_progress=False):   
    
    """Auto-added docstring. See manuscript methods for details."""
    mat_asym=1/(mat_tracer)
    mat_asym[mat_asym==inf]=0
    mat_asym=np.nan_to_num(mat_asym, 0)
    
    mat_sym=1/symmetrize_tracer_mat(mat_tracer)
    mat_sym[mat_sym==inf]=0
    mat_sym=np.nan_to_num(mat_sym, 0)
    
    tracer_only_shuffle=[]
    dmri_only_shuffle=[]
    both_shuffle=[]
    
    for i in range(n_trials):
        mat_shuffled_dmri=shuffle_tracer_preserve_row_idx(mat_dmri)
        if print_progress and i%100==0:
            print('trial ', i)
        
        mat_shuffled_tracer=shuffle_tracer_preserve_row_idx(mat_tracer)
        mat_shuffled_tracer=1/mat_shuffled_tracer
        mat_shuffled_tracer[mat_shuffled_tracer==inf]=0
        mat_shuffled_tracer=np.nan_to_num(mat_shuffled_tracer, 0)
        
        all_shuffled=mat_shuffled_dmri*mat_shuffled_tracer            
       
        #
        #
        if nodeid==-1:
            tracer_only_shuffle.append(directed_mspl(nx.from_numpy_array(mat_dmri*mat_shuffled_tracer,
                                                                        create_using=nx.DiGraph)))
            dmri_only_shuffle.append(directed_mspl(nx.from_numpy_array(mat_shuffled_dmri*mat_asym,
                                                                        create_using=nx.DiGraph)))  
            both_shuffle.append(directed_mspl(nx.from_numpy_array(all_shuffled, 
                                                                                 create_using=nx.DiGraph)))
        else:
            tracer_only_shuffle.append(directed_mspl(nx.from_numpy_array(mat_dmri*mat_shuffled_tracer, 
                                                                        create_using=nx.DiGraph), 
                                                            nodeid=nodeid))
            dmri_only_shuffle.append(directed_mspl(nx.from_numpy_array(mat_shuffled_dmri*mat_asym,
                                                                        create_using=nx.DiGraph),
                                                              nodeid=nodeid))  
            both_shuffle.append(directed_mspl(nx.from_numpy_array(all_shuffled, 
                                                                                 create_using=nx.DiGraph), 
                                                             nodeid=nodeid))    
            
    if nodeid==-1:
        return directed_mspl(nx.from_numpy_array(mat_dmri*mat_sym, create_using=nx.DiGraph)),\
                directed_mspl(nx.from_numpy_array(mat_dmri*mat_asym, create_using=nx.DiGraph)),\
                tracer_only_shuffle, dmri_only_shuffle, both_shuffle

    else:
        return directed_mspl(nx.from_numpy_array(mat_dmri*mat_sym, create_using=nx.DiGraph), nodeid=nodeid),\
               directed_mspl(nx.from_numpy_array(mat_dmri*mat_asym, create_using=nx.DiGraph), nodeid=nodeid),\
               tracer_only_shuffle, dmri_only_shuffle, both_shuffle
               
def process_flow(_dir, species, idx, n_shuffles=100):    
    """Auto-added docstring. See manuscript methods for details."""
    mat_dmri, mat_tracer, gm_labels=get_cmats(root_dir=_dir, 
                                              species=species,
                                             exclude_indices=[15, 21, 24,25,26])
    
    #
    # exclude_indices=[15, 21, 24,25,26] 
    #       - Exclude Claustrum (15), Substantia Nigra (21), Subthalamus (24)
    #       - Homologous in principle, but segmentation limitations in practice. 
    #       - See for more details "A common cross-species atlas of cortical gray matter" https://doi.org/10.1101/2025.09.08.675002
    #
    # Also exclude cerebellum (25) and brainstem (26) for cortical/subcortical gray matter focus. 
    #
        
    np.fill_diagonal(mat_dmri, 0)
    np.fill_diagonal(mat_tracer, 0)    
    vols=[]    
    vol_white=pd.read_table(_dir+species+'_wmstats.txt', index_col=0).T['volume (mm^3)'].sum()    
    vol_white=vol_white**(1/3)    
    mat_dmri=mat_dmri/vol_white    
    
    spl_dmri, spl_dmri_tracer, spls_rand_asymm_dmri, spls_rand_shuffled_dmri, spls_rand_allshuffle = \
                            get_avg_spls(mat_dmri,
                                        mat_tracer,                                                                   
                                        n_trials=n_shuffles, 
                                        print_progress=True)

    return spl_dmri, \
        spl_dmri_tracer, \
        spls_rand_asymm_dmri, \
        spls_rand_shuffled_dmri, \
        spls_rand_allshuffle
