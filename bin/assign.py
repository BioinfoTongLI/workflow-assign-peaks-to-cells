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
import os
import pandas as pd
import pickle
import fire
import numpy as np
from dask_image.imread import imread


def main(peaks, cells, stem="out"):
    with open(cells, "rb") as handle:
        cells = pickle.load(handle)
        cell_ids = np.array(list(cells.keys()))
        cells = np.array(list(cells.values()))
    with open(peaks, "rb") as handle:
        trees = pickle.load(handle)
        if len(trees) == 1:
            trees = trees[0]

    spot_counts = {}
    cell_ys, cell_xs, ys, xs, genes, cell_indexes = [], [], [], [], [], []

    for gene in trees:
        indices_inside = trees[gene].query(cells, predicate="contains")
        target_cells = cells[indices_inside[0]]
        target_cell_ids = cell_ids[indices_inside[0]]
        cell_indexes += target_cell_ids.tolist()
        target_spots = trees[gene].geometries.take(indices_inside[1])
        for cell, spot in zip(target_cells, target_spots):
            ys.append(spot.y)
            xs.append(spot.x)
            genes.append(gene)
            cell_ys.append(cell.centroid.y)
            cell_xs.append(cell.centroid.x)

    spots_df = pd.DataFrame(
        {"y": ys, "x": xs, "genes": genes, "ID": cell_indexes, "cell_y": cell_ys, "cell_x": cell_xs}
    ).set_index("ID")
    spots_df.to_csv(f"{stem}_assigned_peaks.csv")

    count_df = spots_df.pivot_table(index='ID', columns='genes', aggfunc='size', fill_value=0)
    count_df.to_csv(f"{stem}_peak_counts.csv")

    cell_coordinates_df = spots_df.groupby("ID")[["cell_y", "cell_x"]].first()
    cell_coordinates_df.to_csv(f"{stem}_cell_centroids.csv")


def naive_assign(label_image: str, peaks:str, stem: str):
    import json
    if os.path.isdir(peaks):
        # peaks is a folder, assuming to be Xenium
        with open(f"{peaks}/experiment.xenium", "r") as json_file:
            md = json.load(json_file)
        label = imread(label_image).compute()
        transcripts = pd.read_csv(f"{peaks}/transcripts.csv.gz")
        ys = (transcripts["y_location"] / md["pixel_size"]).values.astype(int)
        xs = (transcripts["x_location"] / md["pixel_size"]).values.astype(int)
    elif peaks.endswith(".wkt"):
        from shapely import from_wkt
        with open(peaks, "r") as handle:
            peaks = from_wkt(handle.read())
        ys, xs = [(geom.y, geom.x) for geom in peaks.geoms]
    else:
        raise ValueError("Unrecognized file format")
    # This is to be amended to add cell type info
    cell_ids_for_transcripts = label[0, ys, xs]
    transcripts["id_with_expansion"] = cell_ids_for_transcripts
    transcripts.to_csv(f"{stem}_with_new_cell_ids.csv")
    

if __name__ == "__main__":
    options = {
        "STRtree": main,
        "naive": naive_assign,
        "version": "0.0.1",
    }
    fire.Fire(main)
