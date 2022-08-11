#!/usr/bin/env python

import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib

from nilearn import image as img
import meshplot as mp

from nibabel.cifti2.cifti2 import Cifti2Image

logging.basicConfig(format="%(asctime)s [SGACC THRESHOLD]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)

def threshold_dscalar(dscalar, percentile):
    # Load in the input cifti image
    cifti_img = nib.load(dscalar)

    # Get data from cifti image into array
    cifti_data = cifti_img.get_fdata()

    # Get threshold value (25th percentile of negative values)
    threshold_val = np.percentile(cifti_data[cifti_data < 0], percentile)

    # Copy the data into new array for thresholding
    threshold_data = cifti_data.copy()

    # Threshold the data copy
    threshold_data[threshold_data > threshold_val] = 0

    # Convert the thresholded data back into an image
    new_img = Cifti2Image(threshold_data, header=cifti_img.header,
                         nifti_header=cifti_img.nifti_header)

    # new_img.to_filename(thresholded_dscalar)

    return new_img

def main():
    parser = argparse.ArgumentParser(description="Given a dscalar image "
                                     "and a percentage, performs "
                                     "percentile-based thresholding")
    
    parser.add_argument('dscalar',
                        type=str,
                        help='Path of dscalar to threshold')

    parser.add_argument('percentile',
                        type=float,
                        help='Float of percentage to threshold dscalar by')

    parser.add_argument('output_file',
                        type=str,
                        help='Path of file to output thresholded dscalar to')

    args = parser.parse_args()
    f_dscalar = args.dscalar
    percentile = args.percentile
    output = args.output_file

    logging.info("Thresholding dscalar...")
    thresholded = threshold_dscalar(f_dscalar, percentile)
    thresholded.to_filename(output)

if __name__ == '__main__':
    main()