#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

params.peaks = [
    [0, "/nfs/team283_imaging/HZ_HLB/playground_Tong/HZ_HLB_hindlimb_20220130_63x_fine_tune/decoded/out_opt_flow_registered_decoded_df.tsv"] //path to decoded pandas.DataFrame
    ]
params.target_col = "Name" //gene name column name
params.separator = "\\t" //separator

params.labels = [
    [0, "/nfs/team283_imaging/HZ_HLB/playground_Tong/HZ_HLB_hindlimb_20220130_63x_fine_tune/HZ_HLB_KR0105_C59-FLEG_Nucleus_b1G_b1A_b1T_b1C_Meas6_A2_F1T1_max.ome_label_expanded.tif"]
    ] // path to label image

params.tilesize_x = 700
params.tilesize_y = 700
params.out_dir = "./test"


process Get_shapely_objects {
    echo true
    cache "lenient"
    conda projectDir + "/conda.yaml"
    /*publishDir params.out_dir, mode:'copy'*/

    input:
    tuple val(id), path(lab)
    val(target_col)
    val(separator)

    output:
    tuple val(stem), path("*_shapely.pickle")

    script:
    stem = id
    """
    label_to_shapely.py -label "${lab}"
    """
}


process Get_grid {
    echo true
    cache "lenient"
    conda projectDir + "/conda.yaml"
    /*publishDir params.out_dir, mode:'copy'*/

    input:
    path(peak)
    val(target_col)
    val(separator)
    val(tilesize_x)
    val(tilesize_y)

    output:
    tuple val(stem), path("*_shapely.pickle")

    script:
    stem = peak.baseName
    """
    generate_grid.py --stem "${stem}" --tilesize_x ${tilesize_x} --tilesize_y ${tilesize_y} --csv_in "${peak}" --target_ch "${target_col}" --sep "${separator}"
    """
}


process Build_STR_trees_per_channel {
    echo true
    conda projectDir + "/conda.yaml"
    /*storeDir params.out_dir*/
    /*publishDir params.out_dir, mode:"copy"*/

    input:
    tuple val(id), path(peak)
    val(target_col)
    val(separator)

    output:
    tuple val(id), path("str_peaks.pickle")

    script:
    """
    str_indexing.py -peak ${peak} -target_ch "${target_col}" -sep "${separator}"
    """
}


process Assign {
    /*echo true*/
    conda projectDir + "/conda.yaml"
    /*publishDir params.out_dir, mode:'copy'*/
    storeDir params.out_dir + "/peaks_in_cells"

    input:
    tuple val(stem), path(cells), path(peaks)

    output:
    path("${stem}_assigned_peaks.csv"), emit: peaks_in_cells
    path("${stem}_summary.csv"), emit: peaks_in_cells_summary
    tuple val(stem), path("${stem}_peak_counts.csv"), emit: peaks_counts
    path("${stem}_cell_centroids.csv"), emit: cell_centroids

    script:
    """
    assign.py -cells "$cells" -peaks "$peaks" -stem "$stem"
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


process Shapely_to_label {
    tag "${tiles}"
    echo true

    conda projectDir + "/conda.yaml"
    publishDir params.out_dir, mode:'copy'

    input:
    path(tiles)

    output:
    path("*_label.tif")

    script:
    """
    shapely_to_label.py --tiles_p "$tiles"
    """
}


process to_h5ad {
    tag "${countTable}"
    echo true

    conda projectDir + "/conda.yaml"
    publishDir params.out_dir, mode:'copy'

    input:
    /*path("*_assigned_peaks.csv"), emit: peaks_in_cells*/
    tuple val(stem), path(countTable)
    path(centroids)

    output:
    path("*.h5ad")

    script:
    """
    count_table_2_h5ad.py --countTable ${countTable} --centroids ${centroids} --stem ${stem}
    """
}




workflow {
    Get_shapely_objects(channel.from(params.labels), params.target_col, params.separator)
    _assign(Get_shapely_objects.out)
}


workflow to_grid {
    Get_grid(params.peaks, params.target_col,
        params.separator, params.tilesize_x, params.tilesize_y)
    _assign(Get_grid.out)
}


workflow _assign {
    take: shaply_objs_with_stem
    main:
        Build_STR_trees_per_channel(channel.from(params.peaks),
            params.target_col, params.separator)
        Assign(shaply_objs_with_stem.join(Build_STR_trees_per_channel.out))
        to_h5ad(Assign.out.peaks_counts, Assign.out.cell_centroids)
}
