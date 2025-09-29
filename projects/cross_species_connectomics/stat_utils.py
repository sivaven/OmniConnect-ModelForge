################################################################################
# File: stat_utils.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Statistical utilities. Provides bootstrapping, Kruskal–Wallis testing,
#              permutation testing, and p-value adjustment methods used across analyses.
################################################################################
from matplotlib import ticker
from itertools import combinations
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.stats import kruskal, rankdata
from statsmodels.stats.multitest import multipletests
import numpy as np
import scipy

def bootstrap_ci(data, num_bootstrap=1000, alpha=0.05, rng=None):
    """Auto-added docstring. See manuscript methods for details."""
    rng = np.random.default_rng(rng)
    n = len(data)
    boot_means = [np.mean(rng.choice(data, size=n, replace=True)) for _ in range(num_bootstrap)]
    return np.percentile(boot_means, [100*alpha/2, 100*(1-alpha/2)])
                
# statistics and p-values for pairwise mean differences using bootstrap
def bootstrap_statistic_and_p_value(data1, data2, num_bootstrap=1000, statistic=np.mean):
    """Auto-added docstring. See manuscript methods for details."""
    observed_statistic = statistic(data1) - statistic(data2)
    bootstrap_statistics = []
    combined_data = np.concatenate([data1, data2])  

    for _ in range(num_bootstrap):
        resampled = np.random.choice(combined_data, size=len(combined_data), replace=True)
        resampled_data1 = resampled[:len(data1)]
        resampled_data2 = resampled[len(data1):]
        bootstrap_statistic = statistic(resampled_data1) - statistic(resampled_data2)
        bootstrap_statistics.append(bootstrap_statistic)
    
    bootstrap_statistics = np.array(bootstrap_statistics)
    p_value = np.mean(np.abs(bootstrap_statistics) >= np.abs(observed_statistic))
    
    return observed_statistic, p_value

def cohens_d(a, b):
    """Auto-added docstring. See manuscript methods for details."""
    a = np.asarray(a); b = np.asarray(b)
    m1, m2 = a.mean(), b.mean()
    s1, s2 = a.std(ddof=1), b.std(ddof=1)
    n1, n2 = len(a), len(b)
    # pooled SD
    s_p = np.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2))
    d = (m1 - m2) / s_p
    return d

def hedges_g(a, b):
    # small-sample correction to Cohen's d
    """Auto-added docstring. See manuscript methods for details."""
    n1, n2 = len(a), len(b)
    d = cohens_d(a, b)
    J = 1 - 3/(4*(n1+n2)-9)
    return J * d
    
def dunn_test(groups, alpha=0.05, method='fdr_bh'):
    """Auto-added docstring. See manuscript methods for details."""
    all_data = np.concatenate(groups)
    ranks = rankdata(all_data)
    group_ranks = [ranks[sum(len(g) for g in groups[:i]):sum(len(g) for g in groups[:i+1])] for i in range(len(groups))]
    n = [len(g) for g in groups]
    N = sum(n)
    mean_ranks = [np.mean(r) for r in group_ranks]
    
    comparisons = list(combinations(range(len(groups)), 2))
    z_scores = []
    p_values = []
    max_contributors = []
    for i, j in comparisons:
        diff = mean_ranks[i] - mean_ranks[j]
        se = np.sqrt((N * (N + 1) / 12) * (1 / n[i] + 1 / n[j]))
        z = diff / se
        z_scores.append(z)
        p = 2 * (1 - scipy.stats.norm.cdf(abs(z)))
        p_values.append(p)
    
        # Identify observation contributing most to difference
        rank_diffs = np.abs(group_ranks[i] - group_ranks[j].mean())  # Deviation from mean
        max_contributors.append(groups[i][np.argmax(rank_diffs)])
        
    reject, p_corrected, _, _ = multipletests(p_values, alpha=alpha, method=method)
    return comparisons, z_scores, p_corrected, reject, max_contributors