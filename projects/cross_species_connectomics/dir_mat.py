################################################################################
# File: dir_mat.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh
# Date: 2025-09-27
# Description: Construct Directed Connectivity Matrix Given:
#               - Gray Matter Atlas file (nii.gz)
#               - Tracer Injection and Projection Densities as nii.gz files             
################################################################################

import numpy as np
import pandas as pd
import glob
import nibabel as nb
import sys

def get_p_overlaps(nifti1, val1, nifti2, val2):    
    p = np.sum(nifti2[nifti1==val1]==val2) / (np.sum(nifti2[nifti2==val2]==val2))
    return p

if __name__ == '__main__':   
    sys.stdout.write('reading list of files...')
    seed_files=glob.glob('/aba_mouse/allsites/injection_densities/*injection*')
    end_files=glob.glob('/aba_mouse/allsites/projection_densities/*projection*')
    atlas_file='/aba_mouse/allsites/mouse_CHA.nii.gz'

    nifti_file=nb.load(atlas_file)
    nifti_data=nifti_file.get_fdata()
    unique_num_labels=np.unique(nifti_data)

    grey_labels=pd.read_table(atlas_file.replace('.nii.gz', '.txt'), index_col=0, header=None, delim_whitespace=' ')

    sys.stdout.write('Sorting...')
    seed_files.sort()
    end_files.sort()
    
    sys.stdout.write('Processing...')    
    dir_mat=np.zeros((len(grey_labels), len(grey_labels)))
    
    for f_i in range(len(seed_files)):
        seed_file=seed_files[f_i]
        end_file=end_files[f_i]
        
        if not seed_file.split('/')[-1][0:9]==end_file.split('/')[-1][0:9]:
            sys.stdout.write(seed_file+'\t'+end_file)
        
        #if f_i%10==0:
        #    print('Processing seed,end files #', f_i)
            
        seed_nifti_data=np.round(nb.load(seed_file).get_fdata())
        end_nifti_data=np.round(nb.load(end_file).get_fdata())
        
        
        for i in unique_num_labels[1:]:
            src_p=get_p_overlaps(nifti_data, i, seed_nifti_data, 1.0)
            if src_p>0:
                for j in unique_num_labels[1:]:
                    tar_p=get_p_overlaps(nifti_data, j, end_nifti_data, 1.0)
                    dir_mat[int(i)-1][int(j)-1]+=src_p*tar_p

    df=pd.DataFrame(dir_mat, columns=grey_labels[1].values, index=grey_labels[1].values)
    df.to_csv('/aba_mouse/allsites/mouse_CHA__tracer.csv')
    
