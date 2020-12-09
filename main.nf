#!/usr/bin/env/ nextflow

params.labels = "/nfs/team283_imaging/AC_LNG/playground_Tong/stardist_segs/*.tif"
params.peaks = "/nfs/team283_imaging/AC_LNG/playground_Tong/anchor-spot_location/*.csv"
params.masks = ""
params.outdir = "/nfs/team283_imaging/AC_LNG/playground_Tong/assigned_spots"

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
    //todo
    label_paths.combine(peak_paths, by:0).combine(mask_paths, by:0)
        .set{to_assign}

}else{
    label_paths
        .combine(peak_paths, by:0)
        .set{to_assign}
}

process assign {
    /*errorStrategy 'ignore'*/
    echo true
    conda 'opencv scikit-image pandas shapely imagecodecs'
    publishDir params.outdir, mode:'copy'

    //maxForks 1

    input:
    tuple stem, lab, peak, mask from to_assign

    output:
    path "*_assigned_peaks.csv" into peaks_in_cells
    path "*_summary.csv" into peaks_in_cells_summary
    path "*_peak_counts.csv" into peaks_counts

    script:
    """
    python ${workflow.projectDir}/assign.py -stem "${stem}" -label "$lab" -peak "$peak" -mask "${mask}"
    """
}

