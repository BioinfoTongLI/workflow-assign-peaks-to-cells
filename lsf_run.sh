#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/lustre/scratch126/cellgen/team283/NXF_WORK/'

TMP_NF_WORK="$MOUNT_POINT/${USER}_assign_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run /lustre/scratch126/cellgen/team283/tl10/workflow-assign-peaks-to-cells/main.nf \
	-params-file $1 \
	-profile lsf \
	-resume
