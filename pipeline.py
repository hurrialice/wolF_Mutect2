
import json
import pandas as pd

from wolf import Task, Workflow, output_helpers

GATK_docker_image = {
  "image": "broadinstitute/gatk",
  "tag": "4.1.4.1",
  "user": "root"
}

# read the input.json file for configuration
class Mutect2(Workflow):

    def workflow(self, 
      t_bam_cloud_path, n_bam_cloud_path,
      pairName):
      
      with open('mutect2.inputs.json', 'r') as f:
        config  = json.load(f)
      
      
      self.split_intervals = Task(
        name = "split_intervals",
        inputs = {
          "ref_fasta": config["Mutect2.ref_fasta"],
          "scatter_count" : config['Mutect2.scatter_count'],
          "wes_intervals" : config['Mutect2.intervals'],
          "user_project": config['user_project']
        },
        script = [
          "set -euxo pipefail",
          """
          export GATK_LOCAL_JAR="/root/gatk.jar"
          """,
          
          """/gatk/gatk SplitIntervals -R ${ref_fasta} -L ${wes_intervals} \
            --scatter-count ${scatter_count} -O intervals"""
        ],
        outputs = {
          "subintervals" : "intervals/*-scattered.interval_list"
        },
        overrides={
          "wes_intervals": None,
          "ref_fasta" : None
        },
        docker = GATK_docker_image
      )
      
      self.get_sample_name = Task(
        name = "get_sample_name",
        inputs = {
          "tumor_bam" : t_bam_cloud_path,
          "normal_bam" : n_bam_cloud_path,
          "user_project" : config["user_project"]
        },
        script = [
          "set -euxo pipefail", 
          """
          export GATK_LOCAL_JAR="/root/gatk.jar"
          """,
          "/gatk/gatk GetSampleName -I ${normal_bam} -O normal_name.txt --gcs-project-for-requester-pays ${user_project}",
          "/gatk/gatk GetSampleName -I ${tumor_bam} -O tumor_name.txt --gcs-project-for-requester-pays ${user_project} "
        ],
        outputs = {
          "tumor_name" : ("tumor_name.txt", output_helpers.read_file),
          "normal_name" : ("normal_name.txt", output_helpers.read_file)
        },
        overrides={
          "tumor_bam": None,
          "normal_bam": None,
          "user_project": None
        },
        docker = GATK_docker_image
      )
      
      
      self.M2 = Task(
        name="scatter_m2",
        inputs={
          "ref_fasta": config['Mutect2.ref_fasta'],
          "gnomad": config['Mutect2.gnomad'],
          "default_pon": config['Mutect2.pon'],
          "variants_for_contamination": config['Mutect2.variants_for_contamination'],
          "tumor_bam": t_bam_cloud_path,
          "normal_bam": n_bam_cloud_path,
          "tumor_name": self.get_sample_name.get_output("tumor_name"),
          "normal_name": self.get_sample_name.get_output("normal_name"),
          "extra_args": config['Mutect2.m2_extra_args'],
          "user_project": config['user_project'],
          "command_mem": config['Mutect2.command_mem'],
          "interval" : self.split_intervals.get_output("subintervals")
        },
        script = [
          "set -euxo pipefail", 
          """
          export GATK_LOCAL_JAR="/root/gatk.jar"
          """,
          """
          /gatk/gatk --java-options "-Xmx${command_mem}" Mutect2 \
            -R ${ref_fasta} \
            -I ${tumor_bam} -tumor ${tumor_name} \
            -I ${normal_bam} -normal ${normal_name} \
            --germline-resource ${gnomad} \
            -pon ${default_pon} \
            -L  ${interval} \
            -O {scatter.vcf} \
            --f1r2-tar-gz f1r2.tar.gz \
            --gcs-project-for-requester-pays ${user_project}\
            ${extra_args}
          """,
          
          """
          /gatk/gatk --java-options "-Xmx${command_mem}" GetPileupSummaries -R ~{ref_fasta}\
            -I ${tumor_bam} --interval-set-rule INTERSECTION -L ${intervals} \
            -V ${variants_for_contamination} -L ~{variants_for_contamination} \
            -O tumor-pileups.table
          """,
          
          """
          /gatk/gatk --java-options "-Xmx${command_mem}" GetPileupSummaries -R ~{ref_fasta}\
            -I ${normal_bam} --interval-set-rule INTERSECTION -L ${intervals} \
            -V ${variants_for_contamination} -L ${variants_for_contamination}\
            -O normal-pileups.table
          """
        ],
        outputs = {
          "scatter-vcf": "scatter.vcf*",
          "tpile": "tumor-pileups.table",
          "npile": "normal-pileups.table",
          "f1r2": "f1r2.tar.gz"
        },
        
        resources = {
            "mem" : "6G",
            "cpus-per-task" : 1
        },
        overrides={
          "ref_fasta": None,
          "variants_for_contamination": None,
          "gnomad": None,
          "normal_bam": None,
          "tumor_bam": None,
          "default_pon": None
        },
        docker = GATK_docker_image,
        dependencies = [self.split_intervals, self.get_sample_name]
      )
      # self.calc_contam = Task(
      #   name = "calculate_contamination",
      #   inputs={
      #     "all_tumor_pile_input": self.M2.get_output("tpile").
      #   }
      # )
      





# END PARAMETERS }}}

#
# workflow definitions {{{
# 

##
## TEMPLATE {{{
##
#
##
## define resources
#TASK_mem = "3"
#
#
#
##
## define task
#TASK_task = Task(
#  name = "",
#  inputs = {
#  },
#  outputs = {
#  },
#  script = [
#  ],
#  resources = {
#    "mem" : str(TASK_mem) + "G"
#  },
#  dependencies = 
#)
#
## }}}