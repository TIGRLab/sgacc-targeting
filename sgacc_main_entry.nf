nextflow.preview.dsl=2

include { getArgumentParser } from "./lib/args"

parser = getArgumentParser(
    title: "SGACC TMS Pipeline",
    description: "An end-to-end pipeline for integrating fMRI data into TMS target\
 derivation and optimization",
    scriptName: "${workflow.scriptName}".toString(),
    note: "Configuration arguments should be defined in a .nf.config file (use -c arg), or\
 a .json file that can be used (use -params-file arg)"
)

parser.addRequired("--fmriprep_output",
    "Path to fmriprep output directory",
    params.fmriprep_output.toString(),
    "FMRIPREP_DIRECTORY")

parser.addRequired("--ciftify_output",
    "Path to ciftify output directory",
    params.ciftify_output.toString(),
    "CIFTIFY_DIRECTORY")

parser.addRequired("--out",
    "Path to output directory",
    params.out.toString(),
    "OUTPUT_DIR")

parser.addOptional("--subjects",
    "Path to subject text file containing 1 BIDS subject/line",
    "SUBJECT_FILE")

parser.addOptional("--num_cpus",
    "Maximum number of threads to use when submitting jobs [Default: $params.num_cpus]",
    "NUM_CPUS")

parser.addOptional("--cache_dir",
    "Create a cache directory to store intermediate results to speed up reruns",
    "CACHE_DIR")

missingArgs = parser.isMissingRequired()
missingConfig = parser.isMissingConfig()

if (params.help) {
    print(parser.makeDoc())
    System.exit(0)
}

if (missingArgs || missingConfig) {
    log.error("Missing required parameters")
    missingArgs.each{ log.error("Missing ${it}") }
    missingConfig.each{ log.error("Missing ${it}") }
    print(parser.makeDoc())
    System.exit(1)
}

include {weightfunc_wf} from "./main.nf" params(params)

log.info("Fmriprep Directory: $params.fmriprep_output")
log.info("Ciftify Directory: $params.ciftify_output")
log.info("Output Directory: $params.out")
if (params.subjects) {
    log.info ("Subject list file provided: $params.subjects")
}


input_channel = Channel.fromPath("$params.fmriprep_output/sub-*", type: 'dir')
                    .map{i -> i.getBaseName()}

if (params.subjects){
    subjects_channel = Channel.fromPath(params.subjects)
                            .splitText(){it.strip()}

    input_channel = input_channel.join(subjects_channel)
}

if (!params.rewrite){
    out_channel = Channel.fromPath("$params.out/sub-*", type: 'dir')
                    .map{o -> [o.getBaseName(), "o"]}
                    .ifEmpty(["", "o"])

    input_channel = input_channel.join(out_channel, remainder: true)
                        .filter{it.last() == null}
                        .map{i,n -> i}
}


process publish_coordinates{

    publishDir path: "${params.out}/sgacc_targeting/${sub}", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(coordinates)

    output:
    tuple val(sub), path(coordinates)

    shell:
    '''
    #!/bin/bash
    echo "Transferring !{coordinates} to boonstim/!{sub} folder..."
    '''
}


process publish_qc{

    publishDir path: "${params.out}/targeting_qc/${sub}", \
               mode: 'move', \
               overwrite: true

    input:
    tuple val(sub),\
    path(target_qc_png), path(target_qc_html)

    output:
    tuple val(sub), path(target_qc_png), path(target_qc_html)

    shell:
    '''
    #!/bin/bash
    echo "Transferring !{target_qc_png} to boonstim/!{sub} folder..."
    '''
}

sgacc_input = input_channel.map { sub -> [
            sub,
            "${params.fmriprep_output}/${sub}",
            "${params.ciftify_output}/${sub}",
        ]}

workflow {
    weightfunc_wf(sgacc_input)

    publish_coordinates(weightfunc_wf.out.coordinate)
    publish_qc(weightfunc_wf.out.target_qc)
}