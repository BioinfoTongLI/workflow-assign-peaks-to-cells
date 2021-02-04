#!/usr/bin/env/ nextflow

params.labels = "/nfs/team283_imaging/AC_LNG/playground_Tong/stardist_segs/*.tif"
params.peaks = "/nfs/team283_imaging/AC_LNG/playground_Tong/anchor-spot_location/*.csv"
params.masks = ""
params.outdir = "/nfs/team283_imaging/AC_LNG/playground_Tong/assigned_spots"
params.target_col = "Channel"
params.separator = ","


Channel.fromPath(params.labels)
    .map{it -> [file(file(it).baseName).baseName, it]}
    .set{label_paths}
Channel.fromPath(params.peaks)
    .map{it -> [file(file(it).baseName).baseName, it]}
    .set{peak_paths}

if (params.masks != ""){
    Channel.fromPath(params.masks)
        .map{it -> [file(file(it).baseName).baseName, it]}
        .set{mask_paths}
}else{
    Channel.fromPath(params.labels)
        .map{it -> [file(file(it).baseName).baseName, ""]}
        .set{mask_paths}
}

label_paths
    .combine(mask_paths, by:0)
    .set{to_assign}

process Label_image_to_shapely {
    echo true
    conda baseDir + "/conda.yaml"
    publishDir params.outdir, mode:'copy'

    input:
    tuple stem, lab, mask from to_assign

    output:
    tuple val(stem), file("cell_shapely.pickle") into cell_shapely

    script:
    """
    python ${baseDir}/label_to_shapely.py -label "$lab" -mask "${mask}"
    """
}


process Build_STR_trees_per_channel {
    echo true
    conda baseDir + "/conda.yaml"

    input:
    tuple stem, peak from peak_paths

    output:
    tuple val(stem), file("str_peaks.pickle") into str_peaks

    script:
    """
    python ${baseDir}/str_indexing.py -peak ${peak} -target_ch "${params.target_col}" -sep "${params.separator}"
    """
}


process Assign {
    echo true
    conda baseDir + "/conda.yaml"
    publishDir params.outdir, mode:'copy'
    /*storeDir params.outdir*/

    //maxForks 1

    input:
    tuple stem, cells, peaks from cell_shapely.combine(str_peaks, by:0)

    output:
    file "*_assigned_peaks.csv" into peaks_in_cells
    file "*_summary.csv" into peaks_in_cells_summary
    file "*_peak_counts.csv" into peaks_counts
    file "*_cell_centroids.csv" into cell_centroids

    script:
    """
    python ${baseDir}/assign.py -stem "${stem}" -cells "$cells" -peaks "$peaks"
    """
}


process Cell_filtering {
    echo true
    conda "/home/ubuntu/.conda/envs/assign-peaks-to-cells"
    publishDir params.outdir, mode:'copy'

    input:
    file assigned_peaks from peaks_in_cells
    file peak_counts_in_cell from peaks_counts
    file centroids from cell_centroids

    output:
    file "${stem_assigned_peaks}_thresholded.tsv"
    file "${stem_peak_counts}_thresholded.tsv"
    file "${stem_centroid}_thresholded.tsv"

    script:
    stem_assigned_peaks = file(assigned_peaks).baseName
    stem_peak_counts = file(peak_counts_in_cell).baseName
    stem_centroid = file(centroids).baseName
    """
    python ${baseDir}/cell_filtering.py -assigned_peaks "${assigned_peaks}" -peak_counts_in_cells ${peak_counts_in_cell} -centroids ${centroids} -threshold_n_spots 15 -assigned_peaks_stem ${stem_assigned_peaks} -peak_counts_stem ${stem_peak_counts} -centroid_stem ${stem_centroid}
    """
}

