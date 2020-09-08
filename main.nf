#!/usr/bin/env/ nextflow

mag = '40x'
params.img_path = "/nfs/RV_END/playground_Tong/Full_organoid_segs/assembled_" + mag + "/*.tif"
params.peak_folders_path = "/nfs/RV_END/playground_Tong/Full_organoid_segs/" + mag + "_peaks/"
params.mask_folder_path = "/nfs/RV_END/playground_Tong/Full_organoid_segs/Filtered_" + mag + "/"
img_paths = Channel.fromPath(params.img_path)

process assign {
    errorStrategy 'ignore'
    echo true
    conda 'opencv scikit-image pandas shapely imagecodecs'
    publishDir "/nfs/RV_END/playground_Tong/" + mag + "_peaks_in_organoid_cells", mode:'copy'

    //maxForks 1

    input:
	path img_p from img_paths

    output:
	path "*_peaks.csv" into peaks_in_cells
	path "*_summary.csv" into peaks_in_cells_summary
	path "*_counts.csv" into peaks_counts
    """
	python ${workflow.projectDir}/assign.py -img_in $img_p -spot_csv_dir $params.peak_folders_path -mask_dir ${params.mask_folder_path}
    """
}

