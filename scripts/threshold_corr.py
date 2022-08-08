#!/usr/bin/env python
import argparse
import logging
from collections import namedtuple

import numpy as np
import nibabel as nib

from nilearn import image as img
import meshplot as mp

logging.basicConfig(format="%(asctime)s [SGACC THRESHOLD]:  %(message)s",
                    datefmt="%Y-%m-%d %I:%M:%S %p",
                    level=logging.INFO)

def threshold_dscalar(dscalar, percentile):

    cifti = img.load_img(dscalar)

    thresholded_dscalar = img.threshold_img(cifti, threshold)

    return thresholded_dscalar

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
    percentle = args.percentile
    output = args.output_file

    logging.info("Thresholding dscalar...")
    thresholded = threshold_dscalar(f_dscalar, percentile)
    thresholded.to_filename(output)

if __name__ == '__main__':
    main()