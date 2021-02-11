#!/usr/bin/env/ nextflow

params.labels = "/nfs/team283_imaging/JSP_HSS/playground_Tong/gmm_decoding_human_brain_158_DAPI_registered_40x/cell_seg/out.label_expanded.tif"
params.peaks = "/nfs/team283_imaging/JSP_HSS/playground_Tong/gmm_decoding_human_brain_158_DAPI_registered_40x/decoded/out.decoded_df.tsv"
params.masks = ""
params.target_col = "Name"
params.separator = "\t"
params.to_grid = true
params.tilesize_x = 700
params.tilesize_y = 700
params.outdir = "/nfs/team283_imaging/JSP_HSS/playground_Tong/Grid_count_human_brain_158/peaks_assigned_to_grid_" + params.tilesize_x + "_" + params.tilesize_y



Channel.fromPath(params.labels)
    .map{it -> [file(file(it).baseName).baseName, it]}
    .set{label_paths}
Channel.fromPath(params.peaks)
    .map{it -> [file(file(it).baseName).baseName, it]}
    .into{peak_paths; peaks_to_get_grid}

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

process Get_shapely_objects {
    echo true
    conda baseDir + "/conda.yaml"
    publishDir params.outdir, mode:'copy'

    input:
    tuple stem, lab, mask from to_assign
    tuple stem, peak from peaks_to_get_grid
    val tilesize_x from params.tilesize_x
    val tilesize_y from params.tilesize_y

    output:
    tuple val(stem), file("*_shapely.pickle") into cell_shapely

    script:
    if (!params.to_grid){
        """
        python ${baseDir}/label_to_shapely.py -label "$lab" -mask "${mask}"
        """
    } else {
        """
        python ${baseDir}/generate_grid.py -stem "${stem}" -tilesize_x ${tilesize_x} -tilesize_y ${tilesize_y} -csv_in "${peak}" -target_ch "${params.target_col}" -sep "${params.separator}"
        """
    }
}


process Build_STR_trees_per_channel {
    echo true
    conda baseDir + "/conda.yaml"
    /*storeDir params.outdir*/
    publishDir params.outdir, mode:"copy"

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

