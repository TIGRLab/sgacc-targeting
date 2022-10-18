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


def main():
    parser = argparse.ArgumentParser(description="Generate a QC visualization"
                                " of the selected coordinate on the left hemisphere")
    parser.add_argument("left_gii", help="Left hemisphere gifti", type=str)
    parser.add_argument("dscalar", help="Correlation map", type=str)
    parser.add_argument("coordinate", help="Text file of coordinates", type=str)
    parser.add_argument("qc_img", help="QC image to output")

    args = parser.parse_args()

    gifti = nib.load(args.left_gii)
    l_verts, l_trigs = gifti.agg_data(('pointset', 'triangle'))

    coords = np.loadtxt(args.coordinate, delimiter=' ')

    cifti_img = nib.load(args.dscalar)
    cifti_data = cifti_img.get_fdata()[0]
    cifti_left = cifti_data[:32492]

    face = np.zeros((l_trigs.shape[0], l_trigs.shape[1] + 1), dtype=int)
    face[:, 0] = 3
    face[:,1:] = l_trigs
    face = face.flatten()

    surf = pv.PolyData(l_verts, face)
    surf.point_data['curvature'] = cifti_left

    pv.start_xvfb()

    p = pv.Plotter(polygon_smoothing=True, off_screen=True)
    p.add_mesh(surf)

    sphere = pv.Sphere(center=coords, radius=3)
    p.add_mesh(sphere, color="red")
    p.enable_anti_aliasing()

    p.camera_position = "xz"
    p.camera.azimuth = 260
    p.camera.elevation = 10
    
    p.screenshot(args.qc_img)
    p.export_html(args.qc_img + ".html")


if __name__ == '__main__':
    main()
