import sys
sys.path.append(".")
import pipeline

from wolf import *

#
# Example code for running the WES pipeline on a single pair

with pipeline.WES_workflow(
	# use NFS image preloaded with reference files
	conf = {
		"nfs_image" : "cga-wes-pipeline",
		"nfs_image_project" : "broad-getzlab-workflows"
	}
) as w:
	w.run(
	# paths to tumor/normal BAMs and BAIs in cloud storage
	t_bam_cloud_path = "gs://...",
	n_bam_cloud_path = "gs://...",
	t_bai_cloud_path = "gs://...",
	n_bai_cloud_path = "gs://...",

	# a string for the name of the pair under analysis (used for naming output files)
	pairName = "pair_name",

	# strings for the names of the tumor/normal sample under analysis (used for naming output files)
	caseName = "tumor_sample_name",
	ctrlName = "normal_sample_name"
	)

#
# To run on multiple pairs, simply issue multiple run commands

# For example, if your pairset is in the Pandas dataframe pairset_dataframe:

with pipeline.WES_workflow(
  conf = {
    "nfs_image" : "cga-wes-pipeline",
    "nfs_image_project" : "broad-getzlab-workflows"
  }
) as w:
	for pair in pairset_dataframe.iterrows():
		w.run(
			t_bam_cloud_path = pair["t_bam_cloud_path"],
			n_bam_cloud_path = pair["n_bam_cloud_path"],
			t_bai_cloud_path = pair["t_bai_cloud_path"],
			n_bai_cloud_path = pair["n_bai_cloud_path"],

			pairName = pair["pair_name"],

			caseName = pair["tumor_sample_name"],
			ctrlName = pair["normal_sample_name"]
		)
