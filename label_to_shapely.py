#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import argparse
from shapely.geometry import Polygon
import cv2
import numpy as np
from scipy import ndimage
import tifffile as tf
import pickle


def get_shapely(label):
    """
    Borrowed from Cellpose
    get outlines of masks as a list to loop over for plotting
    """
    polygons = {}
    slices = ndimage.find_objects(label)
    for i, bbox in enumerate(slices):
        if not bbox:
            continue
        cur_cell_label = i + 1
        msk = label[bbox[0], bbox[1]] == cur_cell_label
        cnts, _ = cv2.findContours(
            msk.astype(np.uint8), cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE
        )
        # Multiple split objects is possible, find the largest one
        largest_index = 0
        if cnts[largest_index].shape[0] <= 3:
            continue  # Cell roi too small (only 2-pixel)
        pt = Polygon((cnts[largest_index] + [bbox[1].start, bbox[0].start]).squeeze())
        # print(pt)
        polygons[cur_cell_label] = pt
    return polygons


def main(args):
    label = tf.imread(args.label)
    if args.mask:
        mask = tf.imread(args.mask)
        masked_label = label * (mask == 2)
    else:
        print("Mask doesn't exist, skipped")
        masked_label = label
    with open("cell_shapely.pickle", "wb") as handle:
        pickle.dump(get_shapely(masked_label), handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-label", type=str, required=True)
    parser.add_argument("-mask", type=str)

    args = parser.parse_args()

    main(args)
