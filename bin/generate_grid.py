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
import fire
import numpy as np
from shapely.strtree import STRtree
from shapely.geometry import Polygon, Point
import dask.dataframe as dd
import dask
import pickle


def generate_tiles(min_x, max_x, min_y, max_y, tilesize_x, tilesize_y):
    grid_x = np.arange(min_x, max_x, tilesize_x)
    grid_x_start = grid_x[:-1]
    grid_x_end = grid_x[1:] + 1
    grid_y = np.arange(min_y, max_y, tilesize_y)
    grid_y_start = grid_y[:-1]
    grid_y_end = grid_y[1:] + 1
    count = 1
    shapely_cells = {}
    for x_ind in range(len(grid_x_start)):
        for y_ind in range(len(grid_y_start)):
            pt = Polygon(
                [
                    (grid_x_start[x_ind], grid_y_start[y_ind]),
                    (grid_x_end[x_ind], grid_y_start[y_ind]),
                    (grid_x_end[x_ind], grid_y_end[y_ind]),
                    (grid_x_start[x_ind], grid_y_end[y_ind]),
                ]
            )
            shapely_cells[count] = pt
            count += 1
    return shapely_cells


def main(stem, csv_in, target_ch, sep, tilesize_x, tilesize_y, x_col="x_int", y_col="y_int"):
    df = dd.read_csv(csv_in, sep=sep, dtype={'Code': 'object'})
    df.columns = map(str.lower, df.columns)
    target_ch = target_ch.lower()
    if target_ch != "":
        mask = (
            (df.loc[:, target_ch] != "background")
            & (df.loc[:, target_ch] != "infeasible")
            & (~df.loc[:, target_ch].isna())
        )
        assigned_df = df[mask]
    else:
        assigned_df = df
    max_x, min_x, max_y, min_y = dask.compute(
        assigned_df[x_col].max() + tilesize_x,
        assigned_df[x_col].min() - tilesize_x,
        assigned_df[y_col].max() + tilesize_y,
        assigned_df[y_col].min() - tilesize_y,
    )
    with open("%s_shapely.pickle" % stem, "wb") as handle:
        pickle.dump(
            generate_tiles(min_x, max_x, min_y, max_y, tilesize_x, tilesize_y),
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )


if __name__ == "__main__":
    fire.Fire(main)
