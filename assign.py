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
import pickle


def main(args):
    with open(args.cells, "rb") as handle:
        cells = pickle.load(handle)

    with open(args.peaks, "rb") as handle:
        trees = pickle.load(handle)

    spot_counts = {}
    cell_centroids = {}
    ys, xs, chs, cell_indexes = [], [], [], []
    for cell_index in cells:
        cell = cells[cell_index]
        cell_centroids[cell_index] = {"y": cell.centroid.y, "x": cell.centroid.x}

        current_counts = {}
        for ch in trees:
            potential_inside = trees[ch].query(cell)
            true_in = [sp for sp in potential_inside if cell.contains(sp)]
            for sp in true_in:
                ys.append(sp.y)
                xs.append(sp.x)
                chs.append(ch)
                cell_indexes.append(cell_index)
            current_counts[ch] = len(true_in)
        spot_counts[cell_index] = current_counts

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
    parser.add_argument("-peaks", type=str, required=True)
    parser.add_argument("-cells", type=str, required=True)

    args = parser.parse_args()

    main(args)
