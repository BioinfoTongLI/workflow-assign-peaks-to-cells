#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Assign peaks in .csv to labelled image
"""
import argparse
import pandas as pd
import tifffile as tf
import numpy as np
import cv2

# from skimage.segmentation import find_boundaries
import scipy
from shapely.strtree import STRtree
from shapely.geometry import Polygon, Point
from scipy.ndimage import labeled_comprehension
import re
import os
from pathlib import Path


def get_shapely(label):
    """
    Borrowed from Cellpose
    get outlines of masks as a list to loop over for plotting
    """
    polygons = {}
    slices = scipy.ndimage.find_objects(label)
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


def get_STRtree_per_channel(spots_df, ch_col_name):
    trees = {}
    spots_df.columns = spots_df.columns.str.lower()
    if ch_col_name != "":
        for i, spots in spots_df.groupby(ch_col_name):
            points = [
                Point(spots.loc[ind, "x"], spots.loc[ind, "y"]) for ind in spots.index
            ]
            tree = STRtree(points)
            trees[i] = tree
    else:
        points = [
            Point(spots_df.loc[ind, "x"], spots_df.loc[ind, "y"])
            for ind in spots_df.index
        ]
        trees["anchor"] = STRtree(points)
    return trees


def main(args):
    label = tf.imread(args.label)
    if args.mask:
        mask = tf.imread(args.mask)
        masked_label = label * (mask == 2)
    else:
        print("Mask doesn't exist, skipped")
        masked_label = label
    # print(len(np.unique(label, return_counts=True)[0]))
    # print(len(np.unique(masked_label, return_counts=True)[0]))
    cells = get_shapely(masked_label)

    # Get the STRtree of spots
    spots_df = pd.read_csv(args.peak, sep=args.sep)
    trees = get_STRtree_per_channel(spots_df, args.target_ch.lower())

    spot_counts = {}
    cell_centroids = {}
    ys, xs, chs, cell_indexes = [], [], [], []
    for cell_index in cells:
        cell = cells[cell_index]
        current_counts = {}
        spots_in = {}
        for ch in trees:
            potential_inside = trees[ch].query(cell)
            true_in = [sp for sp in potential_inside if cell.contains(sp)]
            spots_in[ch] = [(sp.y, sp.x) for sp in true_in]
            for sp in true_in:
                ys.append(sp.y)
                xs.append(sp.x)
                chs.append(ch)
                cell_indexes.append(cell_index)
            current_counts[ch] = len(true_in)
        spot_counts[cell_index] = current_counts
        cell_centroids[cell_index] = {"y": cell.centroid.y, "x": cell.centroid.x}
    spots_df = pd.DataFrame(
        {"y": ys, "x": xs, "ch": chs, "ID": cell_indexes}
    ).set_index("ID")
    spots_df.to_csv("%s_assigned_peaks.csv" % args.stem)
    count_df = pd.DataFrame(spot_counts).T
    count_df.to_csv("%s_peak_counts.csv" % args.stem)
    centroid_df = pd.DataFrame(cell_centroids).T
    centroid_df.to_csv("%s_cell_centroids.csv" % args.stem)
    n_total = [trees[i]._n_geoms for i in trees]
    summary = pd.DataFrame(
        {
            "Sum": count_df.sum(),
            "Total_per_ch": n_total,
            "Percentage": count_df.sum() / n_total,
        }
    )
    summary.to_csv("%s_summary.csv" % args.stem)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-stem", type=str, required=True)
    parser.add_argument("-label", type=str, required=True)
    parser.add_argument("-peak", type=str, required=False)
    parser.add_argument("-mask", type=str)
    parser.add_argument("-target_ch", type=str, default="Channel")
    parser.add_argument("-sep", type=str, default="\t")

    args = parser.parse_args()

    main(args)
