#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""

"""
import fire
import shapely
import pickle
import numpy as np
from tifffile import imwrite


def main(tiles_p):
    with open(tiles_p, "rb") as handle:
        tiles = pickle.load(handle)
    last_tile_ind = list(tiles.keys())[-1]

    x_min, y_min, x_max, y_max = tiles[last_tile_ind].bounds
    label = np.zeros((int(y_max), int(x_max)), dtype=np.uint16)
    for tile_index in tiles:
        x_min, y_min, x_max, y_max = tiles[tile_index].bounds
        label[int(y_min):int(y_max), int(x_min):int(x_max)] = tile_index
    imwrite("test_label.tif", label)

if __name__ == "__main__":
    fire.Fire(main)
