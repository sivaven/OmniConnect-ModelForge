################################################################################
# File: vis_utils.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Visualization utilities. Defines plotting functions for connectome evolution,
#              efficiency comparisons, and customized figure aesthetics.
################################################################################

from scipy.interpolate import make_interp_spline
from scipy.stats import sem, t
import scipy.stats
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

def vis_pairwise_stats(stat4vis, 
                       species_names=['Mouse', 'Marmoset', 'Rhesus', 'Human'],
                      xlims=[1,4],
                      ylims=[0,3],
                      save_fig=False,
                       plt_title='',
                      fname_suffix='\\FIGURE4\\fig4_d2.svg',
                       vmin_max=2.5):
    
    """Auto-added docstring. See manuscript methods for details."""
    fig=plt.figure(figsize=(2.25,2.25))
    ax=sns.heatmap(stat4vis, annot=True, fmt=".2f", annot_kws={"fontsize": 8},
                cmap="coolwarm", 
                    vmin=-vmin_max, vmax=+vmin_max,
                xticklabels=species_names, 
                yticklabels=species_names)
    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=7)
    #plt.gca().set_yticklabels(labels, rotation=0, fontsize=10)
    plt.title(plt_title, fontsize=8)
    #plt.xlabel("Species")
    #plt.ylabel("Species")
    plt.xlim([1,4])
    plt.ylim([0,3])

    plt.tight_layout()
    plt.show()
    
def vis_pairwise_stats2(stat4vis, 
                       species_names=['Mouse', 'Marmoset', 'Rhesus', 'Human'],
                      xlims=[1,4],
                      ylims=[0,3],
                      save_fig=False,
                       plt_title='',
                      fname_suffix='',
                       vmin_max=5):
    
    """Auto-added docstring. See manuscript methods for details."""
    fig = plt.figure(figsize=(1.1, 1.3))  
    gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1], hspace=0.3)

    ax = fig.add_subplot(gs[0])
    heatmap = sns.heatmap(stat4vis, annot=True, fmt=".1f", annot_kws={"fontsize": 8},
                          cmap="coolwarm", 
                          vmin=-vmin_max, vmax=+vmin_max,
                          xticklabels=['' for _ in species_names], 
                          yticklabels=['' for _ in species_names],
                          cbar=False, ax=ax)

    # Add horizontal colorbar below
    cbar_ax = fig.add_subplot(gs[1])
    colorbar = fig.colorbar(ax.collections[0], cax=cbar_ax, orientation='horizontal')
    colorbar.ax.tick_params(labelsize=7)

    ax.set_title(plt_title, fontsize=8)
    ax.set_xlim([1, 4])
    ax.set_ylim([0, 3])

    plt.show()

def draw_mean_ci(ax, x, mean_val, lo, hi, color='darkred', lw=.7, ms=3, cap=3, z=6):
    """Auto-added docstring. See manuscript methods for details."""
    lower = max(0.0, mean_val - lo)
    upper = max(0.0, hi - mean_val)
    ax.errorbar(x, mean_val, yerr=[[lower],[upper]],
                fmt='', color=color, linewidth=lw, capsize=cap, zorder=z)
                
def add_break_marks(ax_main, ax_inset, d=0.008, color='k'):
    # top of main axis
    """Auto-added docstring. See manuscript methods for details."""
    km = dict(transform=ax_main.transAxes, color=color, clip_on=False)
    ax_main.plot((-d, +d), (1-d, 1+d), **km)
    ax_main.plot((1-d, 1+d), (1-d, 1+d), **km)
    # bottom of inset axis
    ki = dict(transform=ax_inset.transAxes, color=color, clip_on=False)
    ax_inset.plot((-d, +d), (-d, +d), **ki)
    ax_inset.plot((1-d, 1+d), (-d, +d), **ki)
    
def plot_ev(_all, adj_p, st, s_or_t='', gm_labels=None, species_names=None):
    """Auto-added docstring. See manuscript methods for details."""
    if gm_labels is None:
        raise ValueError("gm_labels must be provided")
    if species_names is None:
        species_names = ['Mouse','Marmoset','Rhesus','Human']
        
    for idx, r in enumerate(gm_labels):    
        _avg = np.array(_all[idx])
        print(r, s_or_t, '\t', adj_p[idx], '\t', st[idx])  

        if adj_p[idx] < 0.05:
            plt.figure(figsize=(1.5, 1.5))

            # Raw traces
            for row in _avg:
                plt.plot(row, '.-', lw=0.5, alpha=.5)

            # Mean and CI
            x = np.arange(len(_avg[0]))
            y = np.mean(_avg, axis=0)
            se = sem(_avg, axis=0)
            ci_range = se * t.ppf(0.975, df=_avg.shape[0] - 1)
            y_lower = y - ci_range
            y_upper = y + ci_range

            # Smoothing
            x_smooth = np.linspace(x.min(), x.max(), 300)
            y_smooth = make_interp_spline(x, y)(x_smooth)
            y_lower_smooth = make_interp_spline(x, y_lower)(x_smooth)
            y_upper_smooth = make_interp_spline(x, y_upper)(x_smooth)

            # Plot mean and CI
            plt.plot(x_smooth, y_smooth, 'k-', lw=1.5)
            plt.fill_between(x_smooth, y_lower_smooth, y_upper_smooth, color='gray', alpha=0.6)

            # Legend and formatting
            plt.plot([], [], lw=0.0001, label=r + '\n(' + s_or_t + ')')
            plt.xticks([0,1,2,3], rotation=90)
            plt.gca().set_xticklabels(['' for _ in species_names], fontsize=9)
            plt.tick_params(axis='y', labelsize=9)
            plt.legend(loc='best', fontsize=8, frameon=False)
            plt.grid(axis='y')

            plt.show()
            
def plot_cost_dict(cost_dicts):    
    """Auto-added docstring. See manuscript methods for details."""
    keys = list(cost_dicts[0].keys())
    cum_values=[0 for k in keys]
    keys.reverse()
    w=.75    
    fig=plt.figure(figsize=(1.9, 3.75))   
    colors=['teal', 'mediumseagreen', 'tab:olive', 'tab:blue']
    for cost_dict, c in zip(cost_dicts, colors):
        values=[cost_dict[k] for k in keys]
        
        key_indices=[len(keys)-list(cost_dict.keys()).index(k) for k in keys]
        
        ax=plt.barh(keys, 
                values, 
                left=cum_values,
                height=w,
                color=c,
                alpha=0.75)
        cum_values=np.add(cum_values, values)

        for i, bars in enumerate(ax.patches[:20]):
            left, bottom, width, height = bars.get_bbox().bounds
            plt.text(left + width / 2, bottom + height / 2, str(key_indices[i]), ha='center', va='center', fontsize=7)
        
    plt.yticks(np.arange(len(keys)), rotation=0)
    
    arrow_labels = [key.replace('__', ' → ') for key in keys]

    plt.gca().set_yticklabels(arrow_labels, fontsize=8)
    plt.gca().set_xticklabels([0,2,4], fontsize=8)     
    plt.xlim([0,4.1]) 
    plt.ylim([-.5,19.5])   
    plt.grid(axis='x', linestyle='--', alpha=.75)
    plt.gca().invert_yaxis()  
    plt.xlabel('Minimum path Cost', fontsize=9)
    #plt.legend(species_names)
    plt.show()
    
def plot_cost_dict_aba_scale(cost_dicts):    
    """Auto-added docstring. See manuscript methods for details."""
    keys = list(cost_dicts[0].keys())
    cum_values=[0 for k in keys]
    keys.reverse()
    w=.75    
    fig=plt.figure(figsize=(3.5,1.25))  
    colors=['black', 'darkred', 'tab:olive', 'tab:blue']
    alphas=[1, .6]
    for cost_dict, c, a in zip(cost_dicts, colors, alphas):
        values=[cost_dict[k] for k in keys]
        
        key_indices=[len(keys)-list(cost_dict.keys()).index(k) for k in keys]
        
        ax=plt.bar(keys, 
                values, 
                bottom=cum_values,
                width=w,
                facecolor='none', edgecolor='black',
                alpha=0.25)
        cum_values=np.add(cum_values, values)       
        
        for i, bars in enumerate(ax.patches[:25]):
            left, bottom, width, height = bars.get_bbox().bounds
            plt.text(left + width / 2, bottom + height / 2, 
                     str(key_indices[i]), ha='center', va='center', fontsize=7, rotation=90, color=c, alpha=a)
        
    plt.xticks(np.arange(len(keys)), rotation=90)
    arrow_labels = [key.replace('__', ' → ') for key in keys]
    plt.gca().set_xticklabels(arrow_labels, fontsize=8)
    plt.gca().set_yticklabels([0,.5, 1], fontsize=8)
    plt.title('Minimum path cost', fontsize=10)
    #plt.xlim([0,2.1]) 
    plt.ylim([0,1.0]) 
    plt.xlim([-.5,24.5])   
    plt.grid(axis='y', linestyle='--', alpha=.75)
    plt.show()