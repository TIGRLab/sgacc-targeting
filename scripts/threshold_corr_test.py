#!/usr/bin/env python

import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib

from nilearn import image as img

from nibabel.cifti2.cifti2 import Cifti2Image

cifti_img = np.array([-10, -8, -7, 0, 1,2,3,4,5])
threshold = 25
print(cifti_img)
print(cifti_img.shape)

dscalar='/projects/smansour/TMS2022/sub-CMP002/tms_pipeline_test/sgacc_targeting_scratch/seed_corr/sub-CMP002_desc-corr.dscalar.nii'
cifti_img = nib.load(dscalar)
print(type(cifti_img))
cifti_data = cifti_img.get_fdata()
wind = 2000
np.set_printoptions(threshold=15)
print(cifti_img)
print(cifti_img.shape)
print(cifti_data[0][wind:wind+10])
print(len(cifti_data[0]))

threshold_val = np.percentile(cifti_data[cifti_data < 0], 25)

cifti_data[cifti_data > threshold_val] = threshold_val
print(cifti_data[0][wind:wind+10])
print(cifti_data.max())
print(cifti_data.min())
print("np")
# print(np.where(cifti_img > threshold_val))
print(threshold_val)
# print(np.where(cifti_data[0] > threshold_val, threshold_val, cifti_data[0]))
print(cifti_data.max())
print(cifti_data.min())
# print(cifti_data[0])

# print(cifti_data)


thresholded_dscalar='/projects/smansour/TMS2022/sub-CMP002/tms_pipeline_test/sgacc_targeting_scratch/threshold_dscalar.nii'

# nib.save(cifti_data, img)

# new_img = nib.Nifti1Image(cifti_data, cifti_img.affine, cifti_img.header)
new_img = Cifti2Image(cifti_data, header=cifti_img.header,
                         nifti_header=cifti_img.nifti_header)

new_img.to_filename(thresholded_dscalar)
print(new_img)

# output_img = nib.load(thresholded_dscalar)
# output_data = output_img.get_fdata()
# print(output_img.shape)