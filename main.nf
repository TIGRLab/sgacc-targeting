nextflow.preview.dsl=2

// To load in a module and forward config parameters to module
// include X from "./PATH/TO/MODULE.nf" params(params)
include {apply_mask as mask_corr} from '../modules/utils.nf' params(params)
include {resample2native_wf as resampleweightfunc_wf} from '../modules/resample2native.nf' params(params)

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

    label 'connectome'
    
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

    label 'bin'

    input:
    tuple val(subject), path(dscalar), val(threshold), 

    output:
    tuple val(subject), path("${subject}_desc-corr_threshold.dscalar.nii"), emit: thresholded

    shell:
    '''
    /scripts/threshold_corr.py !{dscalar} !{threshold} !{subject}_desc-corr_threshold.dscalar.nii
    '''
}

process find_clusters {

    /*
    Generates a cluster report csv based on the input corr map

    Arguments:
        subject (str): Subject ID
        dscalar (Path): Input map to generate clusters from

    Outputs:
        correlation (channel): (subject, corr_map: Path) sgacc correlation map
    */

    label 'connectome'
    
    input:
    tuple val(subject), path(dscalar)

    output:
    tuple val(subject), path("${subject}_desc-clusters.dscalar.nii"), emit: correlation

    shell:
    """
    ciftify_statclust_report ${dscalar} \
        --outputname ${subject}_desc-clusters.dscalar.nii
    """
}

workflow {

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
        seed_corr_input = derivatives.map { sub, fmriprep, ciftify -> [
            sub,
            "${ciftify}/MNINonLinear/" +
            "Results/ses-01_task-rest_desc-preproc/" +
            "ses-01_task-rest_desc-preproc_Atlas_s0.dtseries.nii",
            params.sgacc_template
        ]}
        seed_corr(seed_corr_input)

        // View output
        seed_corr.out.correlation | view

        // 2. thresholding on dscalar
        threshold_corr_input = seed_corr.out.correlation.map {sub, dscalar -> [sub, dscalar, params.threshold]}
        threshold_corr(threshold_corr_input)

        // 3. mask_corr
        mask_corr_input = threshold_corr_input.out.correlation.map {sub, dscalar -> [sub, dscalar, params.dlpfc_mask]}
        mask_corr(mask_corr_input)

        // 4. resample
        resample_inputs = derivatives
                            .spread ( ['L','R'] )
                            .map{subject,folder,h ->   [
                                                subject,
                                                h,
                                                "${ciftify}/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.${h}.sphere.32k_fs_LR.surf.gii" 
                                            ]
                                }
        resampleweightfunc_wf(resample_inputs)

        // 5. clustering
        mask_corr_input = seed_corr.out.correlation.map {sub, dscalar -> [sub, dscalar, params.dlpfc_mask]}
        mask_corr(mask_corr_input)

        // // corr map, Ciftify MSM sphere 
        // seed_corr_mask_input = seed_corr.out.correlation.map {sub, dscalar -> [sub, dscalar, params.dlpfc_mask]}

        // mask_corr(seed_corr_mask_input)
        // weight_func_input = seed_corr.out.correlation.map {sub, dscalar -> [sub, dscalar, params.dlpfc_mask]}

        // resampleweightfunc_wf(mask_corr.out.mask_corr, registration_wf.out.msm_sphere)
        // mapped_inputs = mock_input.spread(['L', 'R']).view().map { subject, ab, h -> [ subject, "${ciftify}/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.${h}.sphere.32k_fs_LR.surf.gii" ]}.view()

        // resample_inputs = derivatives
        //                     .spread ( ['L','R'] )
        //                     .map{subject,folder,h ->   [
        //                                         subject,
        //                                         h,
        //                                         "${ciftify}/${subject}/MNINonLinear/fsaverage_LR32k/${subject}.${h}.sphere.32k_fs_LR.surf.gii" 
        //                                     ]
        //                         }
        // resampleweightfunc_wf(resample_inputs)


    emit:
        // A process which emits the single target coordinate should be here
        coordinate = <PROCESS_NAME>.out.coordinate
}
