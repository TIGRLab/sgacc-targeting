# sgACC Targeting Workflow

![image](https://user-images.githubusercontent.com/54225067/208187504-855cac63-40b5-4dfb-9c10-8d962d99b269.png)

**About**

This is our sgACC targeting pipeline. It's based off of the cluster targeting algorithm described in [this paper from Cash et al.](https://doi.org/10.1002/hbm.25330) Using data in surface (specifically cifti) space, this workflow generates correlation maps seeded from the sgACC, thresholds out the most negative voxels, generates clusters based on the thresholded data, then finally gets the center of mass of the largest cluster. The pipeline also generates visual QC images and html pages for each target to visualize the target on top of the cluster, as shown abovfe.

**Setup**

To use this workflow, you must have nextflow installed on your system. Nextflow is a workflow management service that allows you to make scalable and reproducible workflows for pipelines. It's helpful for our sake because we can write simple workflows via nextflow and run it on our system, and nextflow itself manages all the parallelization. 

You can get started with running it on their [get started page.](https://www.nextflow.io/docs/latest/getstarted.html) You just need nextflow to be able to run on the system, then nextflow handles the pipelines running on the system, and it relies on containers to run each step with software.

On top of this, you'll need containers for the following software in the workflow:
 - ciftify
 - python
 - connectome workbench

These containers will then need to be entered as parameters for the pipeline. These can either be entered as arguments to the nextflow command you end up using, but to avoid bloating that command you can enter it in a json and pass it in as a params file. We have provided an example params.json you can use as a template to get you started. For each parameter in the json, you must enter a value so nextflow can refer to the correct container/config file at runtime.


**Usage**

You can run the sgacc targeting workflow as follows:

```
nextflow sgacc_main_entry.nf \
--ciftify <ciftify_output> --subjects <sub_list> \
-c main.nf.config -c processes.nf.config -params-file params.json \
--out <output_folder> --cache_dir <cache_folder> --sgacc_workflow_directory <sgacc_repo>
```

Here's a quick rundown of the parameters:

```
REQUIRED
	--out	Path to output directory
	--ciftify_output	Path to ciftify output directory
	-c	Path to nextflow config file(s), can use multiple times

OPTIONAL
	--cache_dir	Create a cache directory to store intermediate results to speed up reruns
	--sgacc_workflow_directory	Path to working sgacc workflow directory repo
	--subjects	Path to subject text file containing 1 BIDS subject/line
```

**Outputs**

Once you run the pipeline, you will get an output that looks like the following:
```
output_folder
└── coordinates
    └── sub-{id}
        └── sub-{id}_coordinates.txt
```
Ie, the space separated coordinates for each subject in their own folders.

In addition to this output folder, there is also a cache directory generated for each subject, which saves each output from each intermediate step in the pipeline. This is also where the qc page and image currently lives, ie in this format: 

```
cache_folder
└── visualize_target
    ├── sub-{id}_coordinates.png
    └── sub-{id}_coordinates.html
```
