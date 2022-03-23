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
import pandas as pd
import scanpy as sc


def main(countTable, centroids, stem):
    adata = pd.read_csv(countTable, index_col=0)
    adata = adata.loc[:, ~adata.columns.isin(['background', 'infeasible'])]
    coord = pd.read_csv(centroids, index_col=0)
    adata = sc.AnnData(adata)
    adata.obsm['spatial'] = coord[['x', 'y']].values
    adata.obs['sample'] = stem

    # Remove cells with no mRNA
    adata.obs['total_counts'] = adata.X.sum(1)
    adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(1)

    adata.var['total_counts'] = adata.X.sum(0)
    adata.var['n_cells_by_counts'] = (adata.X > 0).sum(0)
    adata.raw = adata
    adata.write_h5ad(f"{stem}.h5ad")


if __name__ == "__main__":
    fire.Fire(main)
