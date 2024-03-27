#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Get the STRtree of spots
"""
import argparse
import pandas as pd
from shapely.geometry import Point
from shapely.strtree import STRtree
import pickle
import fire


def get_STRtree_per_channel(spots_df, x_col="x_int", y_col="y_int", target_col=None):
    trees = {}
    spots_df.columns = spots_df.columns.str.lower()
    if target_col:
        for i, spots in spots_df.groupby(target_col):
            points = [
                Point(spots.loc[ind, x_col], spots.loc[ind, y_col])
                for ind in spots.index
            ]
            trees[i] = STRtree(points)
    else:
        points = [
            Point(spots_df.loc[ind, x_col], spots_df.loc[ind, y_col])
            for ind in spots_df.index
        ]
        trees["anchor"] = STRtree(points)
    return trees


def main(stem, peak, sep, x_col, y_col, target_col=None):
    spots_df = pd.read_csv(peak, sep=str(sep))
    if target_col in spots_df.columns:
        spots_df = spots_df[
            (spots_df[target_col] != "background")
            & (spots_df[target_col] != "infeasible")
            & (~spots_df[target_col].isna())
            & (spots_df[target_col] != "1.0")
        ]
        str_tree = get_STRtree_per_channel(spots_df, x_col, y_col, target_col.lower()),
    else:
        str_tree = get_STRtree_per_channel(spots_df, x_col, y_col),
    with open(f"{stem}_str_peaks.pickle", "wb") as handle:
        pickle.dump(
            str_tree,
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


if __name__ == "__main__":
    fire.Fire(main)
