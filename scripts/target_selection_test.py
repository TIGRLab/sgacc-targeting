#!/usr/bin/env python

import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib

from nilearn import image as img
# import meshplot as mp

# Dscalar 
dscalar='/projects/smansour/TMS2022/sub-CMP002/tms_pipeline_test/sgacc_targeting_scratch/find_clusters/sub-CMP002_desc-clusters.dscalar.nii'


# Load in the dscalar
cifti_img = nib.load(dscalar)

# Get data from cifti image into array
cifti_data = cifti_img.get_fdata()
print(cifti_data)

# 