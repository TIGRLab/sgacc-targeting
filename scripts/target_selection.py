#!/usr/bin/env python

import argparse
import logging
from collections import namedtuple

import sys

import numpy as np
import nibabel as nib

from nilearn import image as img
import pandas as pd

from nibabel.cifti2.cifti2 import Cifti2Image

CIFTI_LENGTH=32492

logging.basicConfig(format="%(asctime)s [SGACC TARGETING]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)
def load_fdata(dscalar_path):
    return nib.load(dscalar_path).get_fdata()

# def target_selection(dscalar, left_surface, right_surface, sulcal_depth, output_file):
def target_selection(clusters, sulcal, corr_map, left_surface, percentile):
    coordinates = [0, 0, 0]
    # Need to load masked cluster, masked sulcal, cropped corr map, left gifti
    clusters_data = load_fdata(clusters)
    sulcal_data = load_fdata(sulcal)
    corr_map_data = load_fdata(corr_map)
    cifti_img = nib.load(corr_map)

    # Sort the nonzero clusters from largest to smallest
    unique, counts = np.unique(clusters_data, return_counts=True)
    sort_ind = np.argsort(counts)
    sort_ind_rev = sort_ind[::-1]
    cluster_counts = counts[sort_ind_rev]
    cluster_vals = unique[sort_ind_rev]
    
    # Loop over the clusters
    for curr_cluster in cluster_vals:
        # If the current cluster is 0, skip
        if curr_cluster == 0:
            continue
        # Binarize the cluster
        cluster_bin = np.where(clusters_data == curr_cluster, 1.0, 0)

        # Threshold the sulcal map
        sulcal_threshold_val = np.percentile(sulcal_data[sulcal_data > 0], percentile)
        # sulcal_thresholded = np.where(sulcal_data > 0, 1, 0)
        threshold_data = sulcal_data.copy()
        threshold_data[threshold_data < sulcal_threshold_val] = 0

        # Mask the sulcal with the binarized cluster
        sulcal_cluster = threshold_data*cluster_bin

        # Multiply the masked sulcal with the corr map
        corr_mask_cluster = corr_map_data*sulcal_cluster

        # Get the sulcal and corr_map (FOR QC!)
        sulcal_map_qc = corr_map_data*threshold_data

        # Normalize the corr map and get left surface
        norm_corr_mask = corr_mask_cluster/corr_mask_cluster.sum()
        left_corr = norm_corr_mask[0:,CIFTI_LENGTH]

        # Load in the left hemisphere gifti, and calculate the inner product with corr map
        left_gifti = nib.load(left_surface)
        v = left_gifti.agg_data('pointset')
        coordinates = left_corr @ v

        sulcal_map_qc_img = Cifti2Image(sulcal_map_qc, header=cifti_img.header,
                            nifti_header=cifti_img.nifti_header)

        # Return the first average coordinate
        return(coordinates, sulcal_map_qc_img)

    logging.error("Could not find a coordinate from generated clusters, check find_clusters output")
    raise ValueError

def main():
    parser = argparse.ArgumentParser(description="Calculates a centre-of-mass "
                                     "from the largest cluster in the "
                                     "input dscalar file.")
    
    parser.add_argument('clusters',
                        type=str,
                        help='Path of input clusters dscalar file')

    parser.add_argument('sulcal',
                        type=str,
                        help='Path to sulcal height dscalar file')

    parser.add_argument('corr_map',
                        type=str,
                        help='Path to correlation map dscalar file')

    parser.add_argument('left_surface',
                        type=str,
                        help='Path to left surface gifti file')

    parser.add_argument('percentile',
                        type=float,
                        help='Float of percentage to threshold dscalar by')

    parser.add_argument('output_file',
                        type=str,
                        help='Path of file to output coordinates to')

    parser.add_argument('sulcal_map_qc_path',
                        type=str,
                        help='Path of file to output sulcal map qc to')

    args = parser.parse_args()
    f_clusters = args.clusters
    f_sulcal = args.sulcal
    f_corr_map = args.corr_map
    f_left_surface = args.left_surface
    percentile = args.percentile
    output_file = args.output_file
    sulcal_map_qc_path = args.sulcal_map_qc_path

    logging.info("Selecting target...")
    target_coord,sulcal_map_qc_data = target_selection(f_clusters, f_sulcal, f_corr_map, f_left_surface, percentile)
    sulcal_map_qc_data.to_filename(sulcal_map_qc_path)
    np.savetxt(output_file, target_coord)

if __name__ == '__main__':
    main()