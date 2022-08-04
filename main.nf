#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

/*params.target_col = "Channel_Index" //gene name column name*/
params.target_col = "feature_name" //gene name column name
/*params.separator = "\\t"*/
params.separator = ","

/*params.tsv = "/nfs/team283_imaging/NS_DSP/playground_Tong/20220616_RNA_spot_counting/quantifications/label_and_peaks.tsv"*/
params.tsv = "/nfs/team283_imaging/playground_Tong/Xenium_testdata/xenium_prerelease_jun20_mBrain_replicates/mBrain_ff_rep1/transcript_info.csv.gz"

params.tilesize_x = 700
params.tilesize_y = 700
params.out_dir = "output/"


process Get_shapely_objects {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"

    storeDir params.out_dir + "/spot_assignment" //, mode:'copy'

    memory "60 GB"

    input:
    tuple val(stem), path(lab)

    output:
    tuple val(stem), path("${stem}_cell_shapely.pickle")

    script:
    """
    label_to_shapely.py -label "${lab}" -stem ${stem}
    """
}


process Get_grid {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode:'copy'*/

    input:
    path(peak)
    val(target_col)
    val(separator)
    val(tilesize_x)
    val(tilesize_y)
    tuple val(x_col), val(y_col)

    output:
    tuple val(stem), path("*_shapely.pickle")

    script:
    stem = peak.baseName
    """
    generate_grid.py --stem "${stem}" --tilesize_x ${tilesize_x} --tilesize_y ${tilesize_y} --csv_in "${peak}" --target_ch "${target_col}" --sep "${separator}" --x_col ${x_col} --y_col ${y_col}
    """
}


process Build_STR_trees_per_channel {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    storeDir params.out_dir + "/spot_assignment"
    /*publishDir params.out_dir, mode:"copy"*/

    memory "60 GB"

    input:
    tuple val(stem), path(peak)
    val(target_col)
    val(separator)
    tuple val(qc_feature), val(threshold)
    tuple val(x_col), val(y_col)


    output:
    tuple val(stem), path("${stem}_str_peaks.pickle")

    script:
    """
    str_indexing.py -peak ${peak} -target_ch "${target_col}" -sep "${separator}" -stem ${stem} \
        -qc_feature ${qc_feature} -threshold ${threshold} -x_col ${x_col} -y_col ${y_col}
    """
}


process Assign {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    /*publishDir params.out_dir, mode:'copy'*/
    storeDir params.out_dir + "/spot_assignment"

    /*memory "300 GB"*/

    input:
    tuple val(stem), path(cells), path(peaks)

    output:
    path("${stem}_assigned_peaks.csv"), emit: peaks_in_cells, optional: true
    path("${stem}_summary.csv"), emit: peaks_in_cells_summary, optional: true
    tuple val(stem), path("${stem}_peak_counts.csv"), emit: peaks_counts
    path("${stem}_cell_centroids.csv"), emit: cell_centroids, optional: true

    script:
    """
    assign.py -cells "$cells" -peaks "$peaks" -stem "$stem"
    """
}


/*
    Skipped for now and leave the filering to user to perform
*/
process Cell_filtering {
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
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
    debug true

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
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
    /*tag "${countTable}"*/
    debug true

    /*memory "300 GB"*/

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/lustre/scratch117/cellgen/team283/imaging_sifs/assignment.sif':
        'my-docker:latest'}"
    containerOptions "${workflow.containerEngine == 'singularity' ? '--nv':'--gpus all'}"
    publishDir params.out_dir, mode:'copy'

    input:
    /*path("*_assigned_peaks.csv"), emit: peaks_in_cells*/
    tuple val(stem), path(countTable)
    path(centroids)

    output:
    tuple val(stem), path("${stem}.h5ad")

    script:
    """
    count_table_2_h5ad.py --countTable ${countTable} --centroids ${centroids} --stem ${stem}
    """
}


workflow {
    Channel.fromPath(params.tsv)
        .splitCsv(header:true)
        .multiMap{it ->
            labels: [it.stem, it.label]
            peaks: [it.stem, it.peaks]
        }.set{input_files}
    Get_shapely_objects(input_files.labels)
    Build_STR_trees_per_channel(input_files.peaks,
            params.target_col, params.separator, ["Probability", 0.9], ["x_int", "y_int"])
    _assign(Get_shapely_objects.out.join(Build_STR_trees_per_channel.out))
}

workflow to_grid {
    Get_grid(params.tsv, params.target_col, params.separator,
        params.tilesize_x, params.tilesize_y, ["x_location", "y_location"])
    Build_STR_trees_per_channel([file(params.tsv).baseName, params.tsv],
        params.target_col, params.separator, ["qc", 0], ["x_location", "y_location"])
    _assign(Get_grid.out.join(Build_STR_trees_per_channel.out))
}

workflow _assign {
    take: shaply_objs_with_stem
    main:
        Assign(shaply_objs_with_stem)
        to_h5ad(Assign.out.peaks_counts, Assign.out.cell_centroids)
}
