#!/usr/bin/env/ nextflow

mag = '40x'
params.labels = "/nfs/RV_END/playground_Tong/Full_organoid_segs/assembled_" + mag + "/*.tif"
params.peaks = "/nfs/RV_END/playground_Tong/Full_organoid_segs/" + mag + "_peaks/"
params.masks = "/nfs/RV_END/playground_Tong/Full_organoid_segs/Filtered_" + mag + "/"
params.outdir = ""

label_paths = Channel.fromPath(params.labels)
peak_paths = Channel.fromPath(params.peaks)
mask_paths = Channel.fromPath(params.masks)

process assign {
    /*errorStrategy 'ignore'*/
    echo true
    conda 'opencv scikit-image pandas shapely imagecodecs'
    publishDir params.outdir, mode:'copy'

    //maxForks 1

    input:
	path lab from label_paths
	path peak from peak_paths
	path mask from mask_paths

    output:
    path "*_peaks.csv" into peaks_in_cells
    path "*_summary.csv" into peaks_in_cells_summary
    path "*_counts.csv" into peaks_counts
    """
	python ${workflow.projectDir}/assign.py -label $lab -peak $peak -mask $mask
    """
}

