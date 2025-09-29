################################################################################
# File: voxel_count_stats.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Computes tracer coverage for a given atlas regions
################################################################################
import os
import shutil
import sys

import numpy as np
import pandas as pd
import glob
import nibabel as nb

def get_voxel_counts(nifti1, val1, nifti2, val2):     
    """Auto-added docstring. See manuscript methods for details."""
    return np.sum(nifti2[nifti2==val2]==val2), np.sum(nifti1[nifti2==val2]==val1)
	
seed_files=glob.glob('/aba_mouse/allsites/injection_densities/*injection*')
atlas_file='/aba_mouse/allsites/ABAminor.nii.gz'

nifti_file=nb.load(atlas_file)
nifti_data=nifti_file.get_fdata()
unique_num_labels=np.unique(nifti_data)
gray_labels=pd.read_table(atlas_file.replace('.nii.gz', '.txt'), index_col=0, header=None, delim_whitespace=' ')

p_covered=[0 for i in range(len(gray_labels))]
		
df=pd.DataFrame([roi_voxels, p_covered]).T
df.columns=['total_voxels', 'tracer_covered_voxels']
df['percent_covered']=df['tracer_covered_voxels']/df['total_voxels']
df['roi']=gray_labels[1].values



def main():
    sys.stdout.write('Sorting...')
    seed_files.sort()
    for f_i in range(len(seed_files)):   
        if f_i%100==0:
            print('processing seed', f_i)
        roi_voxels=[]
        seed_file=seed_files[f_i]
        seed_nifti_data=np.round(nb.load(seed_file).get_fdata())
        for i in unique_num_labels[1:]:
            roi_vox, roi_seed_vox=get_voxel_counts(seed_nifti_data, 1.0, nifti_data, i)
            roi_voxels.append(roi_vox)
            p_covered[int(i)-1]+=roi_seed_vox
    df.set_index('roi', inplace=True)		
    df.sort_values(by='percent_covered', ascending=False)
    df.to_csv('/aba_mouse/allsites/voxel_count_stat.csv')


if __name__ == '__main__':
    main()
