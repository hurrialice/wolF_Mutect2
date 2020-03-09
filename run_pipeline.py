import sys
sys.path.append(".")
import pipeline

from wolf import *

#
# Example code for running the WES pipeline on a single pair

with pipeline.Mutect2() as w:
	w.run(
	# paths to tumor/normal BAMs and BAIs in cloud storage
	t_bam_cloud_path = "gs://5aa919de-0aa0-43ec-9ec3-288481102b6d-temp/tcga/ACC/DNA/WXS/BI/ILLUMINA/TCGA_MC3.TCGA-OR-A5J7-10A-01D-A29L-10.bam",
	n_bam_cloud_path = "gs://5aa919de-0aa0-43ec-9ec3-288481102b6d-temp/tcga/ACC/DNA/WXS/BI/ILLUMINA/TCGA_MC3.TCGA-OR-A5J7-01A-11D-A29I-10.bam",

	# a string for the name of the pair under analysis (used for naming output files)
	run_name = "ACC_OR-A5J7",
 	# a project to bill for accessing requester-pay buckets
	user_project = "broad-cga-qing-gdac"
	)

#
# To run on multiple pairs, simply issue multiple run commands

# For example, if your pairset is in the Pandas dataframe pairset_dataframe:

# with pipeline.WES_workflow(
#   conf = {
#     "nfs_image" : "cga-wes-pipeline",
#     "nfs_image_project" : "broad-getzlab-workflows"
#   }
# ) as w:
# 	for pair in pairset_dataframe.iterrows():
# 		w.run(
# 			t_bam_cloud_path = pair["t_bam_cloud_path"],
# 			n_bam_cloud_path = pair["n_bam_cloud_path"],
# 			t_bai_cloud_path = pair["t_bai_cloud_path"],
# 			n_bai_cloud_path = pair["n_bai_cloud_path"],

# 			pairName = pair["pair_name"],

# 			caseName = pair["tumor_sample_name"],
# 			ctrlName = pair["normal_sample_name"]
# 		)
