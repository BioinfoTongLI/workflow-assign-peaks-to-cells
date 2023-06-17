#! /bin/sh
#
# build.sh
# Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the MIT license.
#


docker build -t gitlab-registry.internal.sanger.ac.uk/tl10/workflow-assign-peaks-to-cells .
singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/workflow-assign.sif docker-daemon://gitlab-registry.internal.sanger.ac.uk/tl10/workflow-assign-peaks-to-cells:latest
