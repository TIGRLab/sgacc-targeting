nextflow.preview.dsl=2

// To load in a module and forward config parameters to module
// include X from "./PATH/TO/MODULE.nf" params(params)




process seed_corr {
    
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
            params.sgacc_template
        ]}
        seed_corr(seed_corr_input)

        seed_corr.out.correlation | view


}
