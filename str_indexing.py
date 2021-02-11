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
import pysnooper


@pysnooper.snoop()
def get_STRtree_per_channel(spots_df, ch_col_name):
    trees = {}
    spots_df.columns = spots_df.columns.str.lower()
    if ch_col_name != "":
        for i, spots in spots_df.groupby(ch_col_name):
            points = [
                Point(spots.loc[ind, "x"], spots.loc[ind, "y"]) for ind in spots.index
            ]
            trees[i] = STRtree(points)
    else:
        points = [
            Point(spots_df.loc[ind, "x"], spots_df.loc[ind, "y"])
            for ind in spots_df.index
        ]
        trees["anchor"] = STRtree(points)
    return trees


def main(args):
    spots_df = pd.read_csv(args.peak, sep=args.sep)
    # target_ch = "Name" if args.target_ch == "" else args.target_ch
    # spots_df = spots_df[
        # (spots_df[target_ch] != "background")
        # & (spots_df[target_ch] != "infeasible")
        # & (~spots_df[target_ch].isna())
    # ]
    # print(spots_df.shape)
    with open("str_peaks.pickle", "wb") as handle:
        pickle.dump(
            get_STRtree_per_channel(spots_df, args.target_ch.lower()),
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-peak", type=str, required=True)
    parser.add_argument("-target_ch", type=str, default="")
    parser.add_argument("-sep", type=str, default="\t")

    args = parser.parse_args()

    main(args)
