nextflow.preview.dsl=2

// To load in a module and forward config parameters to module
// include X from "./PATH/TO/MODULE.nf" params(params)


process seed_corr {

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
        
        seed_corr_input = derivatives.map { sub, fmriprep, ciftify -> [
            sub,
            "${ciftify}/MNINonLinear/" +
            "Results/ses-01_task-rest_desc-preproc/" +
            "ses-01_task-rest_desc-preproc_Atlas_s0.dtseries.nii",

            // Refers to variables provided in command-line or any config file
            params.sgacc_template
        ]}
        seed_corr(seed_corr_input)

        // View output
        seed_corr.out.correlation | view

    emit:
        // A process which emits the single target coordinate should be here
        coordinate = <PROCESS_NAME>.out.coordinate
}
