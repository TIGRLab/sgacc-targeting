nextflow.preview.dsl=2

// To load in a module and forward config parameters to module
// include X from "./PATH/TO/MODULE.nf" params(params)
// include {apply_mask as mask_corr} from '${projectDir.getParent()}/modules/utils.nf' params(params)
include {apply_mask as mask_corr} from './modules/utils.nf' params(params)

req_param = ["--bids": "$params.bids",
             "--out": "$params.out"]
req_config_param = [
                    "fmriprep": "$params.fmriprep",
                    "ciftify": "$params.ciftify",
                    "connectome": "$params.connectome",
                    "license": "$params.license",
                    "fmriprep_invocation": "$params.fmriprep_invocation",
                    "fmriprep_anat_invocation": "$params.fmriprep_anat_invocation",
                    "fmriprep_descriptor": "$params.fmriprep_descriptor",
                    "ciftify_invocation": "$params.ciftify_invocation",
                    "ciftify_descriptor": "$params.ciftify_descriptor",
                    "weightworkflow" : "$params.weightworkflow"
                   ]

process seed_corr {

    /*
    Generates a correlation map based on sgacc seed mask.

    Arguments:
        subject (str): Subject ID
        dtseries (Path): Path to cifti space timeseries to seed
        sgacc (Path): Path to sgacc roi file

    Outputs:
        correlation (channel): (subject, corr_map: Path) sgacc correlation map
    */

    label 'ciftify'
    
    input:
    tuple val(subject), path(dtseries), path(sgacc)

    output:
    tuple val(subject), path("${subject}_desc-corr.dscalar.nii"), emit: correlation

    shell:
    """
    ciftify_seed_corr ${dtseries} ${sgacc} \
        --outputname ${subject}_desc-corr.dscalar.nii
    """
}

process threshold_corr{

    /*
    Percentile-based thresholding on dscalar file

    Arguments:
        subject (str): Subject ID
        dscalar (Path): Path to dscalar file to threshold
        threshold (float): Float of percentage to threshold dscalar by

    Outputs:
        thresholded (channel): (subject, threshold_dscalar: Path) Thresholded dscalar file
    */

    label 'fieldopt'
    label 'bin'

    input:
    tuple val(subject), path(dscalar), val(threshold)

    output:
    tuple val(subject), path("${subject}_desc-threshold.dscalar.nii"), emit: thresholded

    shell:
    '''
    python /scripts/threshold_corr.py !{dscalar} !{threshold} !{subject}_desc-threshold.dscalar.nii
    '''
}

process find_clusters {

    /*
    Generates a cifti file with nonzero integers for all brainordinates 
    within a large enough cluster, and zeros elsewhere.

    Arguments:
        subject (str): Subject ID
        dscalar (Path): Input map to generate clusters from 
        surface_value_threshold (val): Threshold for surface data values (default: -2.85)
        surface_minimum_area (val): Threshold for surface cluster area, in mm^2 (default: 20)
        volume_value_threshold (val): Threshold for volume data values (default: -2.85)
        volume_minimum_size (val): Threshold for volume cluster size, in mm^3 (default: 20)

    Outputs:
        clusters (channel): (subject, dscalar: Path) sgacc clusters map
    */

    label 'connectome'
    
    input:
    tuple val(subject), path(dscalar),\
    val(surface_value_threshold), val(surface_minimum_area), val(volume_value_threshold), val(volume_minimum_size),\
    path(left_surface), path(right_surface)


    output:
    tuple val(subject), path("${subject}_desc-clusters.dscalar.nii"), emit: clusters

    shell:
    '''
    wb_command -cifti-find-clusters \
        !{dscalar} \
        !{surface_value_threshold} \
        !{surface_minimum_area} \
        !{volume_value_threshold} \
        !{volume_minimum_size} \
        COLUMN \
        !{subject}_desc-clusters.dscalar.nii \
        -left-surface !{left_surface} \
        -right-surface !{right_surface}
    '''
}

process target_selection{

    /*
    Calculate a centre-of-mass from the largest cluster in the input dscalar.

    Arguments:
        subject (str): Subject ID
        dscalar (Path): Path of input clusters dscalar file
        left_surface (Path): Path to left surface
        right_surface (Path): Path to right surface
        sulcal_depth (val): Depth of sulcal to threshold by

    Outputs:
        coordinates (channel): (subject, coordinates: Path) Coordinates of centre of mass target
    */

    label 'fieldopt'
    label 'bin'

    input:
    tuple val(subject), path(dscalar), path(left_surface), path(right_surface), val(sulcal_depth)

    output:
    tuple val(subject), path("${subject}_coordinates.txt"), emit: coordinates

    shell:
    '''
    python /scripts/target_selection.py \
        !{dscalar} \
        !{left_surface} \
        !{right_surface} \
        !{sulcal_depth} \
        !{subject}_coordinates.txt
    '''
}

workflow sgacc_targeting {
    /*
    *   Derivatives tuple (subject: value, fmriprep: path, ciftify: path)
    *       subject: Subject string
    *       fmriprep: Path to subject fMRIPrep folder
    *       ciftify: Path to subject ciftify folder
    */

    take:
        derivatives

    main:
        // 1. seed_corr
        // derivatives | view
        seed_corr_input = derivatives.map { sub, fmriprep, ciftify -> 
                            [
                                sub,
                                "${ciftify}/MNINonLinear/" +
                                "Results/ses-01_task-rest_run-1_desc-preproc/" +
                                "ses-01_task-rest_run-1_desc-preproc_Atlas_s0.dtseries.nii",
                                params.sgacc_template
                            ]}
        seed_corr_input | view
        seed_corr(seed_corr_input)

        // View output
        seed_corr.out.correlation | view

        // 2. thresholding on dscalar
        threshold_corr_input = seed_corr.out.correlation.map {sub, dscalar -> [sub, dscalar, params.threshold]}
        threshold_corr(threshold_corr_input)

        // 3. mask_corr
        mask_corr_input = threshold_corr.out.thresholded.map {sub, dscalar -> [sub, dscalar, params.dlpfc_mask]}
        mask_corr(mask_corr_input)

        // 4. clustering
        cluster_input = mask_corr.out.masked
                                .join(derivatives, by:0)
                                .map { sub, dscalar, fmriprep, ciftify -> 
                                [
                                    sub,
                                    dscalar,
                                    "${params.surface_value_threshold}",
                                    "${params.surface_minimum_area}",
                                    "${params.volume_value_threshold}",
                                    "${params.volume_minimum_size}",
                                    "${ciftify}/T1w/fsaverage_LR32k/" +
                                    "${sub}.L.midthickness.32k_fs_LR.surf.gii",
                                    "${ciftify}/T1w/fsaverage_LR32k/" +
                                    "${sub}.R.midthickness.32k_fs_LR.surf.gii"
                                ]}
        find_clusters(cluster_input)

        // 5. Target Selection
        target_selection_input = seed_corr.out.correlation.map {sub, dscalar -> [sub, dscalar, params.threshold]}
        target_selection_input = find_clusters.out.clusters
                                        .join(derivatives, by:0)
                                        .map { sub, dscalar, fmriprep, ciftify -> 
                                        [
                                            sub,
                                            dscalar,
                                            "${ciftify}/T1w/fsaverage_LR32k/" +
                                            "${sub}.L.midthickness.32k_fs_LR.surf.gii",
                                            "${ciftify}/T1w/fsaverage_LR32k/" +
                                            "${sub}.R.midthickness.32k_fs_LR.surf.gii",
                                            "${params.sulcal_depth}"
                                        ]}
        target_selection(target_selection_input)

    emit:
        // A process which emits the single target coordinate should be here
        coordinate = target_selection.out.coordinates
}
