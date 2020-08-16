#!/usr/bin/env/ nextflow

mag = '40x'
params.img_path = "/nfs/RV_END/playground_Tong/Full_organoid_segs/assembled_" + mag + "/Assembled/*.tif"
params.peak_folders_path = "/nfs/RV_END/playground_Tong/Full_organoid_segs/" + mag + "_peaks/"
img_paths = Channel.fromPath(params.img_path)
// filePairs = Channel.fromFilePairs("/nfs/RV_END/playground_Tong/Full_organoid_segs/**ome_{cyto_diam_90_assembled.tif, peaks.csv}")

process assign {
    errorStrategy 'ignore'
    echo true
    conda 'opencv scikit-image pandas shapely'
    publishDir "/nfs/RV_END/playground_Tong/" + mag + "_peaks_in_organoid_cells", mode:'copy'

    input:
	path img_p from img_paths

    output:
	path "*_assigned_peaks.csv" into peaks_in_cells
    """
	python ${workflow.projectDir}/assign.py -img_in $img_p -spot_csv_dir $params.peak_folders_path
    """
}

