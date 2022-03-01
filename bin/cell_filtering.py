#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Filter cells base on number of assigned spots
"""
import argparse
import pandas as pd
import numpy as np


def main(args):
    assigned_peaks = pd.read_csv(args.assigned_peaks, index_col=0)
    peak_counts = pd.read_csv(args.peak_counts_in_cells, index_col=0)
    centroids = pd.read_csv(args.centroids, index_col=0)

    thresholded_peak_counts = peak_counts[
        peak_counts.sum(axis=1) > args.threshold_n_spots
    ]
    thresholded_peak_counts.to_csv(args.peak_counts_stem + "_thresholded.tsv", sep="\t")
    peaks_mask = [ind in thresholded_peak_counts.index for ind in assigned_peaks.index]
    thresholded_assigned_peaks = assigned_peaks[peaks_mask]
    thresholded_assigned_peaks.to_csv(
        args.assigned_peaks_stem + "_thresholded.tsv", sep="\t"
    )

    centroids_mask = [ind in thresholded_peak_counts.index for ind in centroids.index]
    thresholded_centroids = centroids[centroids_mask]
    thresholded_centroids.to_csv(args.centroid_stem + "_thresholded.tsv", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-assigned_peaks", type=str, required=True)
    parser.add_argument("-peak_counts_in_cells", type=str, required=True)
    parser.add_argument("-centroids", type=str, required=True)
    parser.add_argument("-centroid_stem", type=str, required=True)
    parser.add_argument("-assigned_peaks_stem", type=str, required=True)
    parser.add_argument("-peak_counts_stem", type=str, required=True)
    parser.add_argument("-threshold_n_spots", type=int, default=15)

    args = parser.parse_args()

    main(args)
