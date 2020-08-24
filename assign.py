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


def get_shapely(label):
    """
    Borrowed from Cellpose
    get outlines of masks as a list to loop over for plotting
    """
    polygons={}
    slices = scipy.ndimage.find_objects(label)
    for i, bbox in enumerate(slices):
        if not bbox:
            continue
        cur_cell_label = i + 1
        msk = label[bbox[0], bbox[1]] == cur_cell_label
        _, cnts, _ = \
            cv2.findContours(
                    msk.astype(np.uint8),
                    cv2.RETR_LIST,
                    cv2.CHAIN_APPROX_SIMPLE
                    )
        # Multiple split objects is possible, find the largest one
        largest_index = 0
        if len(cnts) > 1:
            largest_index = np.argmax(
                    [cv2.contourArea(cnt) for cnt in cnts]
                    )
        if cnts[largest_index].shape[0] <= 3:
            continue # Cell roi too small (only 2-pixel)
        pt = Polygon(
                (cnts[largest_index] + [bbox[1].start, bbox[0].start]).squeeze()
                )
        # print(pt)
        polygons[cur_cell_label] = pt
    return polygons


def get_STRtree_per_channel(spots_df):
    trees = {}
    for i, spots in spots_df.groupby('Channel'):
        points = [Point(spots.loc[ind, 'x'],
            spots.loc[ind, 'y'])
                for ind in spots.index]
        tree = STRtree(points)
        trees[i] = tree
    return trees


# def shapely_fast(mark):
    # mark[find_boundaries(mark, mode='inner')] = 0
    # mark[:,[0,-1]], mark[[0,-1],:] = 0, 0   # if we want leave the mark touch the border
    # _, contours, hierarchy = cv2.findContours((mark==0).astype(np.uint8), cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    # print(hierarchy)
    # contours = [contours[i] for i in np.where(hierarchy[0,:,3]>0)[0]] # only inner contours is needed
    # print(contours)


def main(args):
    stem = re.search(r'(.*ome)_.*', args.img_in).group(1)
    label = tf.imread(args.img_in)
    cells = get_shapely(label)

    # Get the STRtree of spots
    spots_df = pd.read_csv("%s/%s_peaks.csv"
            %(args.spot_csv_dir, stem))
    trees = get_STRtree_per_channel(spots_df)

    spot_counts = {}
    ys, xs, chs, cell_indexes = [], [], [], []
    for cell_index in cells:
        cell = cells[cell_index]
        current_counts = []
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
            current_counts.append(len(true_in))
        spot_counts[cell_index] = current_counts
    spots_df = pd.DataFrame({"y":ys, "x":xs, "ch":chs, "ID":cell_indexes}).set_index("ID")
    spots_df.to_csv("%s_peaks.csv" %stem)
    df = pd.DataFrame(spot_counts).T
    df.rename(columns= {0:"ch1", 1:"ch2", 2:"ch3"},
            inplace=True)
    n_total = [trees[i]._n_geoms for i in trees]
    summary = pd.DataFrame({'Sum':df.sum(),
            'Total_per_ch':n_total,
            'Percentage':df.sum()/n_total})
    summary.to_csv("%s_summary.csv" %stem)
    df.to_csv("%s_peak_counts.csv" %stem)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-img_in", type=str,
            required=True)
    parser.add_argument("-spot_csv_dir", type=str,
            required=True)

    args = parser.parse_args()

    main(args)
