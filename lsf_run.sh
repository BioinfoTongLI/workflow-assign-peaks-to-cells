#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

PATH="/software/singularity-v3.6.4/bin/":$PATH
MOUNT_POINT='/lustre/scratch117/cellgen/team283/NXF_WORK/'

DATE_WITH_TIME=`date "+%Y%m%d%H%M"`
TRACE_FILE="$MOUNT_POINT/${USER}_${DATE_WITH_TIME}_trace/registration_trace_${DATE_WITH_TIME}.tsv"
TMP_NF_WORK="$MOUNT_POINT/${USER}_${DATE_WITH_TIME}_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK LSB_DEFAULTGROUP='team283' /lustre/scratch117/cellgen/team283/tl10/nextflow/nextflow -trace nextflow.executor run main.nf \
	-params-file $1 \
	-with-trace $TRACE_FILE \
	-profile lsf \
	#-resume
