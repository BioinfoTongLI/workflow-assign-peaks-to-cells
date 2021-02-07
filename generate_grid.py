#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Assign peaks to grid
"""
import argparse
import numpy as np
from shapely.strtree import STRtree
from shapely.geometry import Polygon, Point
import dask.dataframe as dd
import dask
import pickle
import pysnooper


def generate_cells(min_x, max_x, min_y, max_y):
    grid_x = np.arange(min_x, max_x, args.tilesize_x)
    grid_x_start = grid_x[:-1]
    grid_x_end = grid_x[1:] + 1
    grid_y = np.arange(min_y, max_y, args.tilesize_y)
    grid_y_start = grid_y[:-1]
    grid_y_end = grid_y[1:] + 1
    count = 1
    shapely_cells = {}
    for x_ind in range(len(grid_x_start)):
        for y_ind in range(len(grid_y_start)):
            pt = Polygon([
                (grid_x_start[x_ind], grid_y_start[y_ind]),
                (grid_x_end[x_ind], grid_y_start[y_ind]),
                (grid_x_end[x_ind], grid_y_end[y_ind]),
                (grid_x_start[x_ind], grid_y_end[y_ind]),
                ])
            shapely_cells[count] = pt
            count += 1
    return shapely_cells


def main(args):
    df = dd.read_csv(args.csv_in, sep=args.sep)
    mask = (df[args.target_ch] != "background") & (df[args.target_ch] != "infeasible") & (~df[args.target_ch].isna())
    assigned_df = df[mask]
    max_x, min_x, max_y, min_y = dask.compute(
            assigned_df.X.max() + args.tilesize_x,
            assigned_df.X.min() - args.tilesize_x,
            assigned_df.Y.max() + args.tilesize_y,
            assigned_df.Y.min() - args.tilesize_y
        )
    with open("%s_shapely.pickle" %args.stem, "wb") as handle:
        pickle.dump(generate_cells(min_x, max_x, min_y, max_y), handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-stem", type=str, required=True)
    parser.add_argument("-csv_in", type=str, required=True)
    parser.add_argument("-target_ch", type=str, required=True)
    parser.add_argument("-sep", type=str, required=True)
    parser.add_argument("-tilesize_x", type=int, default=350)
    parser.add_argument("-tilesize_y", type=int, default=350)

    args = parser.parse_args()

    main(args)
