nextflow.preview.dsl=2

process apply_mask {

    label 'connectome'
    input:
    tuple val(sub), path(dscalar), path(mask), val(name_suffix)

    output:
    tuple val(sub), path("${sub}.${name_suffix}.dscalar.nii"), emit: masked

    shell:
    '''
    wb_command -cifti-math \
                "x * (mask > 0)" \
                -var "x" !{dscalar} \
                -var "mask" !{mask} \
                !{sub}.!{name_suffix}.dscalar.nii
    '''
}

process average_corr_maps{

    label 'connectome'

    input:
    tuple val(sub), path("*")

    output:
    tuple val(sub), path("${sub}_merged.dscalar.nii"), emit: merged_dscalar

    shell:
    '''
    find . -mindepth 1 -maxdepth 1 -type l -name "*dscalar.nii" | sort | xargs -I {} \
        echo -cifti {} | xargs \
        wb_command -cifti-average !{sub}_merged.dscalar.nii
    '''

}

process remove_subcortical{

    label 'connectome'

    input:
    tuple val(sub), path(dscalar)

    output:
    tuple val(sub), path("${sub}.correlation_nosubcort.dscalar.nii"), emit: corr_dscalar

    shell:
    '''
    # Split without volume
    wb_command -cifti-separate !{dscalar} \
                COLUMN \
                -metric CORTEX_LEFT L.shape.gii \
                -metric CORTEX_RIGHT R.shape.gii
    # Join
    wb_command -cifti-create-dense-scalar \
                !{sub}.correlation_nosubcort.dscalar.nii \
                -left-metric L.shape.gii \
                -right-metric R.shape.gii
    '''

}