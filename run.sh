#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#

nextflow run main.nf -params-file yamls/FrFr_lung.yaml
nextflow run main.nf -params-file yamls/FrFr_brain.yaml
nextflow run main.nf -params-file yamls/FrFr_breast.yaml
