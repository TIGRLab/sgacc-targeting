#!/usr/bin/env python

import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib

from nilearn import image as img
import meshplot as mp

logging.basicConfig(format="%(asctime)s [SGACC TARGETING]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)

def target_selection(dscalar, left_surface, right_surface, sulcal_depth, output_file):
    coordinates_file = [0, 0, 0]
    # Iteratively,


        # Select the largest cluster

        # Threshold using sulcal depth

        # Calculate the centre-of-mass

        # If the cluster is removed during sulcal depth thresholding, find the next biggest cluster and loop
    
    return coordinates_file

def main():
    parser = argparse.ArgumentParser(description="Calculates a centre-of-mass "
                                     "from the largest cluster in the "
                                     "input dscalar file.")
    
    parser.add_argument('dscalar',
                        type=str,
                        help='Path of input clusters dscalar file')

    parser.add_argument('left_surface',
                        type=str,
                        help='Path to left surface')

    parser.add_argument('right_surface',
                        type=str,
                        help='Path to right surface')

    parser.add_argument('sulcal_depth',
                        type=float,
                        help='Depth of sulcal to threshold by')

    parser.add_argument('output_file',
                        type=str,
                        help='Path of file to output coordinates to')

    args = parser.parse_args()
    f_dscalar = args.dscalar
    f_left_surface = args.left_surface
    f_right_surface = args.right_surface
    sulcal_depth = args.sulcal_depth
    output_file = args.output_file

    logging.info("Selecting target...")
    target_coord = target_selection(f_dscalar, f_left_surface, f_right_surface, sulcal_depth, output_file)
    np.savetxt(output_file, target_coord)

if __name__ == '__main__':
    main()