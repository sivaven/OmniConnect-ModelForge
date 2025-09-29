################################################################################
# File: fig3a.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Generates Figure 3a. Computes bootstrap confidence intervals and visualizations
#              of cross-species path efficiency distributions.
################################################################################

import cross_species_nets as cn
import stat_utils as su
import vis_utils as vu

import numpy as np
import matplotlib.pyplot as plt
from vis_utils import add_break_marks

from matplotlib import ticker
from itertools import combinations

_dir='data/'

dmri_syms=[]
dmri_asyms=[]
dmri_ran_asyms=[]
ran_dmri_asyms=[]
ran_boths=[]

speciess=['mouse',  'marmoset', 'rhesus', 'human']

species_names = ["Mouse", "Marmoset", "Rhesus", "Human"]
pairwise_results = []
pair_ds = {}

def main():
    # collect data per species
    for idx, species in enumerate(speciess):
        dmri_sym, dmri_asym, dmri_ran_asym, ran_dmri_asym, ran_both = cn.process_flow(
            _dir, species, idx, n_shuffles=100 #<- set to 10000
        )
        dmri_syms.append(dmri_sym)
        dmri_asyms.append(dmri_asym)
        dmri_ran_asyms.append(dmri_ran_asym)
        ran_dmri_asyms.append(ran_dmri_asym)
        ran_boths.append(ran_both)

    species_data = [r for r in dmri_asyms]
    null_means_per_species = [np.mean(r, axis=1) for r in dmri_ran_asyms]    

    fig=plt.figure(figsize=(1.75,2))
    gs  = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.1, 1.9], hspace=0.05)
    ax_main  = fig.add_subplot(gs[1])
    ax_inset = fig.add_subplot(gs[0], sharex=ax_main)

    cutoff_main = 20
    pos_real = [1,2,3,4]
    pos_null = [p+0.4 for p in pos_real]
    w_real, w_null = 0.4, 0.2

    bp_r_main = ax_main.boxplot(
        species_data, positions=pos_real, widths=[w_real]*4, patch_artist=True,
        medianprops={"color":"white","linewidth":1},
        flierprops={"marker":'.',"markerfacecolor":'r',"markersize":1}, labels=species_names
    )
    bp_n_main = ax_main.boxplot(
        null_means_per_species, positions=pos_null, widths=[w_null]*4, patch_artist=True,
        medianprops={"color":"white","linewidth":1},
        flierprops={"marker":'.',"markerfacecolor":'r',"markersize":.25}, labels=['']*4
    )

    # inset plots
    bp_r_in = ax_inset.boxplot(
        species_data, positions=pos_real, widths=[w_real]*4, patch_artist=True,
        showmeans=False, medianprops={"color":"white","linewidth":1},
        flierprops={"marker":'.',"markerfacecolor":'r',"markersize":1,  "alpha":0.0},  
        whiskerprops={"alpha":0}, capprops={"alpha":0},
        labels=['']*4
    )
    bp_n_in = ax_inset.boxplot(
        null_means_per_species, positions=pos_null, widths=[w_null]*4, patch_artist=True,
        showmeans=False, medianprops={"color":"white","linewidth":1}, 
        flierprops={"marker":'.',"markerfacecolor":'r',"markersize":.25, "alpha":0.0}, 
        whiskerprops={"alpha":.0}, capprops={"alpha":0}, labels=['']*4
    )

    # styling
    ax_main.yaxis.grid(True, linestyle='--', alpha=0.75)
    ax_inset.yaxis.grid(True, linestyle='--', alpha=0.75)
    ax_main.set_ylim(-2, cutoff_main)
    ax_inset.set_ylim(1e0, 1e2)
    ax_inset.set_yscale("log")
    ax_main.set_yticks([0, 10, 15])
    ax_main.set_yticklabels([0, 10, 15], fontsize=8)    
    ax_main.set_ylabel("Minimum path cost", fontsize=9)
    ax_inset.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=3))
    ax_inset.yaxis.set_minor_locator(ticker.NullLocator())
    ax_inset.set_yticks([1, 10, 100])
    ax_inset.set_yticklabels([ '1', '10', '100'], fontsize=8)

    for b in bp_r_main['boxes']:
        b.set_facecolor('darkseagreen'); b.set_linewidth(.5)
    for b in bp_n_main['boxes']:
        b.set_facecolor('grey'); b.set_edgecolor('black'); b.set_linewidth(.5)
    for b in bp_r_in['boxes']:
        b.set_facecolor('darkseagreen'); b.set_linewidth(.5); b.set_alpha(0.0)
    for b in bp_n_in['boxes']:
        b.set_facecolor('grey'); b.set_linewidth(.5); b.set_alpha(0.0)

    # add mean ± CI on inset
    for i, data in enumerate(species_data, start=1):
        lo, hi = su.bootstrap_ci(np.asarray(data))
        vu.draw_mean_ci(ax_inset, i, float(np.mean(data)), lo, hi, color='darkred')
    for i, null_means in enumerate(null_means_per_species, start=1):
        lo, hi = np.percentile(null_means, [2.5, 97.5])
        vu.draw_mean_ci(ax_inset, i+0.4, float(np.mean(null_means)), lo, hi, color='darkred')

    ax_inset.spines['bottom'].set_visible(False)
    ax_main.spines['top'].set_visible(False)
    ax_inset.tick_params(labelbottom=False)
    ax_main.xaxis.tick_bottom()
    add_break_marks(ax_main, ax_inset)
    ax_main.set_xticks(pos_real)
    ax_main.set_xticklabels([], fontsize=8, rotation=20)
    plt.tight_layout()

    # pairwise stats
    stat4vis=np.zeros((len(species_data), len(species_data)))+np.nan
    for i, data1 in enumerate(species_data):
        for j, data2 in enumerate(species_data):
            if i < j:
                observed_stat, p_value = su.bootstrap_statistic_and_p_value(data1, data2)
                stat4vis[i][j]=observed_stat
                pairwise_results.append((species_names[i], species_names[j], observed_stat, p_value))
                print(f"Statistic for {species_names[i]} vs {species_names[j]}: {observed_stat:.4f}, P-value: {p_value:.4f}")
    for (i,j) in combinations(range(4), 2):
        d = su.cohens_d(species_data[i], species_data[j])
        g = su.hedges_g(species_data[i], species_data[j])    
        pair_ds[(species_names[i], species_names[j])] = dict(d=d, g=g)  
    for k,v in pair_ds.items():
        print(f"{k[0]} vs {k[1]}: Cohen's d = {v['d']:.3f} | Hedges' g = {v['g']:.3f}")
    for i, data1 in enumerate(species_data):
        for j, data2 in enumerate(species_data):
            if i < j:           
                stat4vis[i][j]=pair_ds[(species_names[i], species_names[j])]['d']
    vu.vis_pairwise_stats(stat4vis, save_fig=False, plt_title='', vmin_max=.8)

if __name__ == '__main__':
    main()
