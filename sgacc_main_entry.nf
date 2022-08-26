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

parser.addRequired("--fmriprep",
    "Path to fmriprep output directory",
    params.fmriprep.toString(),
    "FMRIPREP_DIRECTORY")

parser.addRequired("--ciftify",
    "Path to ciftify output directory",
    params.ciftify.toString(),
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

include {sgacc_targeting} from "./main.nf" params(params)

log.info("Fmriprep Directory: $params.fmriprep")
log.info("Ciftify Directory: $params.ciftify")
log.info("Output Directory: $params.out")
if (params.subjects) {
    log.info ("Subject list file provided: $params.subjects")
}


input_channel = Channel.fromPath("$params.fmriprep/sub-*", type: 'dir')
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

    publishDir path: "${params.out}/coordinates/${sub}", \
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

sgacc_input = input_channel.map { sub -> [
            sub,
            "${params.fmriprep}/${sub}",
            "${params.ciftify}/${sub}",
        ]}

workflow {
    sgacc_targeting(sgacc_input)

    publish_coordinates(sgacc_targeting.out.coordinate)
}