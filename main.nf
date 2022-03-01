#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

params.labels = "/nfs/team283_imaging/JSP_HSS/playground_Tong/gmm_decoding_human_brain_158_DAPI_registered_40x/cell_seg/out.label_expanded.tif"
params.peaks = "/nfs/team283_imaging/JSP_HSS/playground_Tong/gmm_decoding_human_brain_158_DAPI_registered_40x/decoded/out.decoded_df.tsv"
params.target_col = "Name"
params.separator = "\\t"

params.to_grid = false
params.tilesize_x = 700
params.tilesize_y = 700
params.out_dir = "./test"


process Get_shapely_objects {
    echo true
    conda projectDir + "/conda.yaml"
    publishDir params.out_dir, mode:'copy'

    input:
    path(lab)
    path(peak)
    val(target_col)
    val(separator)
    val(to_grid)
    val(tilesize_x)
    val(tilesize_y)

    output:
    path("*_shapely.pickle")

    script:
    if (to_grid){
        """
        generate_grid.py -stem "${stem}" -tilesize_x ${tilesize_x} -tilesize_y ${tilesize_y} -csv_in "${peak}" -target_ch "${target_col}" -sep "${separator}"
        """
    } else {
        """
        label_to_shapely.py -label "${lab}"
        """
    }
}


process Build_STR_trees_per_channel {
    echo true
    conda projectDir + "/conda.yaml"
    /*storeDir params.out_dir*/
    publishDir params.out_dir, mode:"copy"

    input:
    path(peak)
    val(target_col)
    val(separator)

    output:
    path("str_peaks.pickle")

    script:
    """
    str_indexing.py -peak ${peak} -target_ch "${target_col}" -sep "${separator}"
    """
}


process Assign {
    /*echo true*/
    conda projectDir + "/conda.yaml"
    publishDir params.out_dir, mode:'copy'

    input:
    path(cells)
    path(peaks)

    output:
    path("*_assigned_peaks.csv"), emit: peaks_in_cells
    path("*_summary.csv"), emit: peaks_in_cells_summary
    path("*_peak_counts.csv"), emit: peaks_counts
    path("*_cell_centroids.csv"), emit: cell_centroids

    script:
    """
    assign.py -cells "$cells" -peaks "$peaks"
    """
}


/*
    Skipped for now and leave the filering to user to perform
*/
process Cell_filtering {
    echo true
    conda projectDir + "/conda.yaml"
    publishDir params.out_dir, mode:'copy'

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
    cell_filtering.py -assigned_peaks "${assigned_peaks}" -peak_counts_in_cells ${peak_counts_in_cell} -centroids ${centroids} -threshold_n_spots 15 -assigned_peaks_stem ${stem_assigned_peaks} -peak_counts_stem ${stem_peak_counts} -centroid_stem ${stem_centroid}
    """
}

workflow {
    Get_shapely_objects(params.labels, params.peaks,
        params.target_col, params.separator, params.to_grid,
        params.tilesize_x, params.tilesize_y)
    Build_STR_trees_per_channel(params.peaks, params.target_col, params.separator)
    Assign(Get_shapely_objects.out, Build_STR_trees_per_channel.out)
}

