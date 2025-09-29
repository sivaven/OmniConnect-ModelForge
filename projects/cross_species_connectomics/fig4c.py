################################################################################
# File: fig4c.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Generates Figure 4c. Processes full path efficiency table to
#              quantify evolutionary changes and visualize pairwise comparison stats.
################################################################################

import cross_species_nets as cn
import stat_utils as su
import vis_utils as vu

import numpy as np
import pandas as pd

sig_rois=['TEM_Inferior', 'TEM_Superior', 'INS_Anterior', 'INS_Anterior', 'OLF_Anterior', 'OLF_Piriform']
s_or_t=['Source', 'Source', 'Source', 'Target', 'Source', 'Target']
fig_fname_sfx='//KW//'

_dir='data/'
df_evolution=pd.read_csv(_dir+"full_path_efficiency_table.csv", index_col=0)

species_names=['Mouse', 'Marmoset', 'Rhesus', 'Human']


def main():
    for roi, src_or_tar in zip(sig_rois, s_or_t):  

        setyticklabels=True
        if roi!='TEM_Superior':
            setyticklabels=False

        stat4vis=np.zeros((len(species_names), len(species_names)))+np.nan

        print(roi, src_or_tar)
        src_avg=[]
        tar_avg=[]        

        for path in df_evolution.index:
            src=path.split('__')[0]
            tar=path.split('__')[1]        
            if src==roi:            
                src_avg.append(df_evolution.loc[path])
            if tar==roi:            
                tar_avg.append(df_evolution.loc[path])


        if src_or_tar=='Source':
            groups=np.array(src_avg).T
            comparisons, z_scores, p_corrected, reject, max_contr = su.dunn_test(groups)
            for (i, j), z, p, r, mc in zip(comparisons, z_scores, p_corrected, reject, max_contr):
                print(f" {species_names[i]} vs {species_names[j]}, \t Z-score: {z:.2f}, \t max-contribution: {mc:.4f}, \tCorrected p-value: {p:.4f}, Reject Null: {r}")
                stat4vis[i][j]=z
            vu.vis_pairwise_stats2(stat4vis,
                               save_fig=True, 
                               plt_title=roi +'\n('+src_or_tar+')',
                               fname_suffix=fig_fname_sfx+'_posthoc_'+roi+'_'+src_or_tar+'.svg')

        if src_or_tar=='Target':
            groups=np.array(tar_avg).T
            comparisons, z_scores, p_corrected, reject, max_contr = su.dunn_test(groups)
            for (i, j), z, p, r, mc in zip(comparisons, z_scores, p_corrected, reject, max_contr):
                print(f" {species_names[i]} vs {species_names[j]}, \t Z-score: {z:.2f}, \t max-contribution: {mc:.4f}, \tCorrected p-value: {p:.4f}, Reject Null: {r}")
                stat4vis[i][j]=z
            vu.vis_pairwise_stats2(stat4vis, 
                               save_fig=True, 
                               plt_title=roi +'\n('+src_or_tar+')',
                               fname_suffix=fig_fname_sfx+'_posthoc_'+roi+'_'+src_or_tar+'.svg')


if __name__ == '__main__':
    main()
