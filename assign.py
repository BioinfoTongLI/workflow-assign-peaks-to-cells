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
from skimage import draw
from skimage.segmentation import find_boundaries
from shapely.strtree import STRtree
from shapely.geometry import Polygon, Point
from scipy.ndimage import labeled_comprehension
import re


def get_shapely(masks):
    """ get outlines of masks as a list to loop over for plotting """
    polygons={}
    for n in np.unique(masks)[1:]:
        mn = masks==n
        _, contours, _ = cv2.findContours(mn.astype(np.uint8), mode=cv2.RETR_LIST, method=cv2.CHAIN_APPROX_SIMPLE)
        cmax = np.argmax([c.shape[0] for c in contours])
        pix = contours[cmax].astype(int).squeeze()
        if len(pix)>4:
            pix=pix[:,::-1]
            pix = draw.polygon_perimeter(pix[:,0], pix[:,1], (mn.shape[0], mn.shape[1]))
            pix = np.array(pix).T[:,::-1]
            polygons[n] = Polygon(pix)
        else:
            polygons[n] = Polygon(np.zeros((0,2)))
    return polygons


def get_STRtree_per_channel(spots_df):
    trees = {}
    for i in spots_df.Channel.unique():
        spots = spots_df[spots_df.Channel == i]
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
    for n_cell in cells:
        cell = cells[n_cell]
        current_counts = []
        for ch in trees:
            potential_inside = trees[ch].query(cell)
            true_in = [sp for sp in potential_inside if cell.contains(sp)]
            current_counts.append(len(true_in))
        spot_counts[n_cell] = current_counts
    df = pd.DataFrame(spot_counts).T
    df.rename(columns= {0:"ch1", 1:"ch2", 2:"ch3"},
            inplace=True)
    df.to_csv("%s_assigned_peaks.csv" %stem)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-img_in", type=str,
            required=True)
    parser.add_argument("-spot_csv_dir", type=str,
            required=True)

    args = parser.parse_args()

    main(args)
