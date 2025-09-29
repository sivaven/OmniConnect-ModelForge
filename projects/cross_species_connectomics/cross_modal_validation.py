################################################################################
# File: cross_modal_validation.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Validation analyses comparing mouse viral tracer data with diffusion MRI tractography,
#              including cross-modal alignment and consistency checks.
################################################################################
import pandas as pd
from scipy import stats
import random
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import glob

#
# Fig. 2c
#
_dir='data/'
trk_op_dir=_dir+'crossmodal_streamline_stats//'
op_sfx='.tt.gz.stat.txt'

stat_files=glob.glob(trk_op_dir+"*"+op_sfx)
inj_ids=[stat_file.split('\\')[-1] for stat_file in stat_files]
inj_ids=[i.split('_')[0] for i in inj_ids]
inj_ids=list(set(inj_ids))

inj_only_seeds=[]
inj_proj_combs=[]
inj_proj_randperm=[]

def main():
    for inj_id in inj_ids:    
        df=pd.read_table(trk_op_dir+inj_id+'_seed'+op_sfx, header=None, index_col=0)
        ntracts=df.loc['number of tracts'][1]
        inj_only_seeds.append(ntracts)    

        df=pd.read_table(trk_op_dir+inj_id+'_seed__'+inj_id+'_end'+op_sfx, header=None, index_col=0)
        ntracts=df.loc['number of tracts'][1]
        inj_proj_combs.append(ntracts)    

        perm_files=glob.glob(trk_op_dir+inj_id+'*')
        proj_perm=inj_id
        for perm_file in perm_files:
            tokens=perm_file.split('\\')[-1].split('_')        
            if len(tokens)>2 and tokens[3]!=tokens[0]:
                proj_perm=tokens[3]
        df=pd.read_table(trk_op_dir+inj_id+'_seed__'+proj_perm+'_end'+op_sfx, header=None, index_col=0)
        ntracts=df.loc['number of tracts'][1]
        inj_proj_randperm.append(ntracts)

    df_cons = pd.DataFrame([inj_only_seeds, inj_proj_combs, inj_proj_randperm]).T
    df_cons.columns = ['seed_only', 'seed_end', 'seed_rand']
    df_cons['frac_end'] = df_cons['seed_end'] / df_cons['seed_only']
    df_cons['frac_rand'] = df_cons['seed_rand'] / df_cons['seed_only']

    # Long-form for violinplot
    df_long = pd.melt(df_cons[['frac_end', 'frac_rand']], var_name='Condition', value_name='Fraction')

    fig, ax = plt.subplots(figsize=(1.25, 1.5))
    vp = sns.violinplot(
        data=df_long,
        x='Condition',
        y='Fraction',
        palette='Set3',
        cut=0,
        inner='box',
        linewidth=1,
        ax=ax
    )

    # Medians based on actual data
    medians = df_long.groupby('Condition')['Fraction'].median()

    # Styling
    for artist in vp.collections:
        artist.set_edgecolor('black')
        artist.set_linewidth(.5)
    for i, (cond, median_val) in enumerate(medians.items()):
        ax.plot(i, median_val, marker='o', markersize=3, color='black', zorder=3)

    ax.text(0.5, 1.09, f'n={len(df_cons)}', ha='center', va='bottom', fontsize=7,
            transform=ax.transAxes)
    ax.set_xticklabels(['seed\n&\nend ', 'seed\n&\nrandom end'], size=8, rotation=0)
    ax.set_ylabel('Fraction of \nSeed-Only Tracts', size=9)
    ax.set_xlabel('')
    ax.set_ylim([-.1, 1.1])
    ax.grid(True, linestyle='-', axis='y', alpha=.75)

    plt.show()

if __name__ == '__main__':
    main()
