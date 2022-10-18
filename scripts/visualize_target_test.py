#!/usr/bin/env python

import argparse
import logging
from collections import namedtuple

import sys

import pyvista as pv
from pathlib import Path
import nibabel as nib
import matplotlib.pyplot as plt
import nilearn.plotting as nplot
import numpy as np
from nibabel.gifti.gifti import GiftiImage


f_left_gii = '/projects/smansour/TMS2022/sub-CMP003/tms_pipeline_test/fmriprep_output/ciftify/sub-CMP003/T1w/fsaverage_LR32k/sub-CMP003.L.midthickness.32k_fs_LR.surf.gii'
dscalar = '/projects/smansour/TMS2022/sub-CMP003/tms_pipeline_test/sgacc_targeting_scratch/crop_corr/sub-CMP003.correlation_nosubcort.dscalar.nii'
f_coordinate = '/KIMEL/tigrlab/projects/smansour/TMS2022/sub-CMP003/tms_pipeline_test/sgacc_targeting_scratch/target_selection/sub-CMP003_coordinates.txt'
qc_img = '/scratch/smansour/scripts/sgacc-dev/sub-CMP003_img'


pv.start_xvfb()

gifti = nib.load(f_left_gii)
l_verts, l_trigs = gifti.agg_data(('pointset', 'triangle'))

coords = np.loadtxt(f_coordinate, delimiter=' ')

cifti_img = nib.load(dscalar)
cifti_data = cifti_img.get_fdata()[0]
cifti_left = cifti_data[:32492]

face = np.zeros((l_trigs.shape[0], l_trigs.shape[1] + 1), dtype=int)
face[:, 0] = 3
face[:,1:] = l_trigs
face = face.flatten()

surf = pv.PolyData(l_verts, face)
surf.point_data['curvature'] = cifti_left


p = pv.Plotter(polygon_smoothing=True, off_screen=True)
p.add_mesh(surf)
# p.start_xvfb()

sphere = pv.Sphere(center=coords, radius=3)
p.add_mesh(sphere, color="red")
p.enable_anti_aliasing()

p.camera_position = "xz"
p.camera.azimuth = 260
p.camera.elevation = 10


p.window_size = [5000, 5000]
p.camera.zoom(0.50)
# p.add_text("Green arrow indicates coil anterior position"
#            "\nScalar map is E-field magnitude (V/m)" + additional_text,
#            color="white",
#            shadow=True,
#            font_size=15)
# p.show()
p.screenshot(qc_img)
p.export_html(qc_img + ".html")