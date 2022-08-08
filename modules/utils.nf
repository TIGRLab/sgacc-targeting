nextflow.preview.dsl=2

process apply_mask {

    label 'connectome'
    input:
    tuple val(sub), path(dscalar), path(mask)

    output:
    tuple val(sub), path("${sub}.masked.dscalar.nii"), emit: masked

    shell:
    '''
    wb_command -cifti-math \
                "x * (mask > 0)" \
                -var "x" !{dscalar} \
                -var "mask" !{mask} \
                !{sub}.masked.dscalar.nii
    '''
}