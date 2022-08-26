nextflow.preview.dsl=2
import groovy.util.FileNameByRegexFinder

// To load in a module and forward config parameters to module
include {apply_mask as mask_corr} from './modules/utils.nf' params(params)
include {apply_mask as mask_sulcal} from './modules/utils.nf' params(params)
include {average_corr_maps as avg_corr_map} from './modules/utils.nf' params(params)
include {average_corr_maps as avg_back_project} from './modules/utils.nf' params(params)
include {remove_subcortical as crop_corr} from './modules/utils.nf' params(params)
include {remove_subcortical as crop_clusters} from './modules/utils.nf' params(params)

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
    tuple val(subject), val(run), path(dtseries), path(sgacc)

    output:
    tuple val(subject), val(run), path("${subject}_${run}_desc-corr.dscalar.nii"), emit: correlation

    shell:
    """
    ciftify_seed_corr ${dtseries} ${sgacc} \
        --outputname ${subject}_${run}_desc-corr.dscalar.nii
    """
}


process seed_corr_back_project {

    /*
    Generates a correlation map based on sgacc seed mask using back projection

    Arguments:
        subject (str): Subject ID
        dtseries (Path): Path to cifti space timeseries to seed
        sgacc (Path): Path to sgacc roi file

    Outputs:
        correlation (channel): (subject, corr_map: Path) sgacc correlation map
    */

    label 'ciftify'
    
    input:
    tuple val(subject), val(run), path(dtseries), path(sgacc)

    output:
    tuple val(subject), val(run), path("${subject}_${run}_desc-corr.dscalar.nii"), emit: correlation

    shell:
    """
    ciftify_seed_corr --weighted ${dtseries} ${sgacc} \
        --outputname ${subject}_${run}_desc-corr.dscalar.nii
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
    val(surface_value_threshold), val(surface_minimum_area),\
    val(volume_value_threshold), val(volume_minimum_size),\
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
        clusters (Path): Path of input clusters dscalar file
        sulcal (Path): Path to sulcal height dscalar file
        corr_map (Path): Path to correlation map dscalar file
        left_surface (Path): Path to left surface gifti file
        sulcal_depth (val): Depth of sulcal to threshold by

    Outputs:
        coordinates (channel): (subject, coordinates: Path) Coordinates of centre of mass target
    */

    label 'fieldopt'
    label 'bin'

    input:
    tuple val(subject), path(clusters), path(sulcal),\
    path(corr_map), path(left_surface), val(sulcal_depth)

    output:
    tuple val(subject), path("${subject}_coordinates.txt"), emit: coordinates

    shell:
    '''
    python /scripts/target_selection.py \
        !{clusters} \
        !{sulcal} \
        !{corr_map} \
        !{left_surface} \
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
        // 1. Seed from the sgacc
        seed_corr_input = derivatives
                            .map{s,f,c ->   [
                                                s,
                                                new FileNameByRegexFinder().getFileNames("${c}",
                                                ".*MNINonLinear/Results/.*(REST|rest).*/.*dtseries.nii")
                                            ]
                                }
                            .transpose()
                            .map{s,run ->   [
                                                s,
                                                ( run =~ /run-[^_]*/ )[0],
                                                run,
                                                params.sgacc_template
                                            ]
                                }
        seed_corr(seed_corr_input)


        // Seed from the dlpfc (back projection method)
        back_project_input = derivatives
                            .map{s,f,c ->   [
                                                s,
                                                new FileNameByRegexFinder().getFileNames("${c}",
                                                ".*MNINonLinear/Results/.*(REST|rest).*/.*dtseries.nii")
                                            ]
                                }
                            .transpose()
                            .map{s,run ->   [
                                                s,
                                                ( run =~ /run-[^_]*/ )[0],
                                                run,
                                                params.sgacc_back_template
                                            ]
                                }
        seed_corr_back_project(back_project_input)

        // 2. Average the correlation maps across runs
        avg_corr_map_inputs = seed_corr.out.correlation.groupTuple( by: 0 , sort: {it}, size: 6 ).map {sub, run, dscalar -> [sub, dscalar]}
        avg_corr_map(avg_corr_map_inputs)

        // Average the correlation maps across runs (back projection method)
        avg_back_project_inputs = seed_corr_back_project.out.correlation.groupTuple( by: 0 , sort: {it}, size: 6 ).map {sub, run, dscalar -> [sub, dscalar]}
        avg_back_project(avg_back_project_inputs)

        // 3. Threshold the correlation map
        threshold_corr_input = avg_back_project.out.merged_dscalar.map {sub, dscalar -> [sub, dscalar, params.threshold]}
        threshold_corr(threshold_corr_input)

        // 4. Crop the thresholded correlation map 
        crop_corr_input = threshold_corr.out.thresholded.map {sub, dscalar -> [sub, dscalar]}
        crop_corr(crop_corr_input)

        // 5. Mask the sulcal with the dlpfc
        mask_sulcal_input = derivatives
                            .map{sub,fmriprep,ciftify -> [sub, "${ciftify}/MNINonLinear/fsaverage_LR32k/${sub}.sulc.32k_fs_LR.dscalar.nii", params.dlpfc_mask, "sulcal_masked"]}
        mask_sulcal(mask_sulcal_input)

        // 6. Generate the clusters from the thresholded correlation map
        cluster_input = crop_corr.out.corr_dscalar
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

        // 7. Mask the clusters with the dlpfc
        mask_corr_input = find_clusters.out.clusters.map {sub, dscalar -> [sub, dscalar, params.dlpfc_mask, "corr_masked"]}
        mask_corr(mask_corr_input)

        // 8. Run the target selection on the masked corr map and masked sulcal
        target_selection_input = mask_corr.out.masked
                                        .join(mask_sulcal.out.masked, by:0)
                                        .join(crop_corr.out.corr_dscalar, by:0)
                                .join(derivatives, by:0)
                                .map { sub, mask_corr, mask_sulcal, crop_corr, fmriprep, ciftify -> 
                                [
                                    sub, mask_corr, mask_sulcal, crop_corr,
                                    "${ciftify}/T1w/fsaverage_LR32k/" +
                                    "${sub}.L.midthickness.32k_fs_LR.surf.gii",
                                    params.sulcal_depth
                                ]}
        target_selection(target_selection_input)

    // Output the target selection coordinate
    emit:
        coordinate = target_selection.out.coordinates
}
