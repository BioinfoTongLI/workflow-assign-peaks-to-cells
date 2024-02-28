#! /bin/sh
#
# build.sh
# Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the MIT license.
#


docker build -t bioinfotongli/workflow-assign-peaks-to-cells .
singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/cache/bioinfotongli-workflow-assign-peaks-to-cells-latest.img docker-daemon://bioinfotongli/workflow-assign-peaks-to-cells:latest
