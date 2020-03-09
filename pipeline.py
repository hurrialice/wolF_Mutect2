
import json
import pandas as pd

from wolf import Task, Workflow, output_helpers

GATK_docker_image = {
  "image": "broadinstitute/gatk",
  "tag": "4.1.4.1",
  "extra_flags" : {"volume" : "/home/qing/.config/gcloud:/etc/gcloud" }
}

# read the input.json file for configuration
class Mutect2(Workflow):

    def workflow(self, 
      t_bam_cloud_path, n_bam_cloud_path,
      user_project):
      
      with open('mutect2.inputs.json', 'r') as f:
        config  = json.load(f)
      
      ref_path = config["ref_local_path"]
      
      
      self.split_intervals = Task(
        name = "split_intervals",
        inputs = {
          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
          "scatter_count" : config['Mutect2.scatter_count'],
          "wes_intervals" : ref_path + config['Mutect2.intervals'],
          "user_project": user_project
        },
        script = [
          "set -euxo pipefail",
          "export CLOUDSDK_CONFIG=/etc/gcloud",
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
          "user_project" : user_project
        },
        script = [
          "set -euxo pipefail",
          "export CLOUDSDK_CONFIG=/etc/gcloud",
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
          "ref_fasta": ref_path + config['Mutect2.ref_fasta'],
          "gnomad": ref_path + config['Mutect2.gnomad'],
          "default_pon": ref_path + config['Mutect2.pon'],
          "variants_for_contamination": ref_path + config['Mutect2.variants_for_contamination'],
          "tumor_bam": t_bam_cloud_path,
          "normal_bam": n_bam_cloud_path,
          "tumor_name": self.get_sample_name.get_output("tumor_name"),
          "normal_name": self.get_sample_name.get_output("normal_name"),
          "extra_args": config['Mutect2.m2_extra_args'],
          "user_project": user_project,
          "command_mem": "4",
          "interval" : self.split_intervals.get_output("subintervals", lambda x: x[0])
        },
        script = [
          "set -euxo pipefail",
          "export CLOUDSDK_CONFIG=/etc/gcloud",
          "echo $interval",
          """
          /gatk/gatk --java-options "-Xmx${command_mem}g" Mutect2 \
            -R ${ref_fasta} \
            -I ${tumor_bam} -tumor ${tumor_name} \
            -I ${normal_bam} -normal ${normal_name} \
            --germline-resource ${gnomad} \
            -pon ${default_pon} \
            -L  ${interval} \
            -O scatter.vcf \
            --f1r2-tar-gz f1r2.tar.gz \
            --gcs-project-for-requester-pays ${user_project}\
            ${extra_args}
          """,
          
          """
          /gatk/gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries -R ${ref_fasta}\
            -I ${tumor_bam} --interval-set-rule INTERSECTION -L ${interval} \
            -V ${variants_for_contamination} -L ${variants_for_contamination} \
            --gcs-project-for-requester-pays ${user_project}\
            -O tumor-pileups.table
          """,
          
          """
          /gatk/gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries -R ${ref_fasta}\
            -I ${normal_bam} --interval-set-rule INTERSECTION -L ${interval} \
            -V ${variants_for_contamination} -L ${variants_for_contamination}\
            --gcs-project-for-requester-pays ${user_project}\
            -O normal-pileups.table
          """
        ],
        outputs = {
          "scatter_vcf": "scatter.vcf*",
          "tpile": "tumor-pileups.table",
          "npile": "normal-pileups.table",
          "f1r2": "f1r2.tar.gz"
        },
        
        resources = {
            "mem" : "6G",
            "cpus-per-task" : 1
        },
        overrides={
          "normal_bam": None,
          "tumor_bam": None
        },
        docker = GATK_docker_image,
        dependencies = [self.split_intervals, self.get_sample_name]
      )
      self.calc_contam = Task(
        name = "calculate_contamination",
        inputs={
          "all_tumor_pile_input": self.M2.get_output("tpile", lambda x: " -I ".join(x)),
          "all_normal_pile_input": self.M2.get_output("npile", lambda x: " -I ".join(x)),
          "seq_dict": ref_path + config['Mutect2.ref_dict']
        },
        script=[
          "set -euxo pipefail",
          """
          /gatk/gatk --java-options -Xmx2g GatherPileupSummaries -I ${all_tumor_pile_input} \
            --sequence-dictionary ${seq_dict} -O tumor_pile.tsv
          """,
          """
          /gatk/gatk --java-options -Xmx2g GatherPileupSummaries -I ${all_normal_pile_input} \
            --sequence-dictionary ${seq_dict} -O normal_pile.tsv
          """,
          """
          /gatk/gatk --java-options -Xmx2g CalculateContamination \
            -I tumor_pile.tsv -O contamination.table \
            --tumor-segmentation segments.table \
            -matched normal_pile.tsv
          """
        ],
        outputs={
          "contamination": "contamination.table",
          "segments": "segments.table",
          "tumor_pile_table" : "tumor_pile.tsv",
          "normal_pile_table" : "normal_pile.tsv"
        },
        docker = GATK_docker_image,
        dependencies=self.M2)
      
      self.merge_M2 = Task(
        name="merge_m2",
        inputs={
          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
          "contamination_table": self.calc_contam.get_output("contamination"),
          "segments_table": self.calc_contam.get_output("segments"),
          "all_f1r2_input": self.M2.get_output("f1r2", lambda x: " -I ".join(x)),
          "all_vcf_input": self.M2.get_output("scatter_vcf", lambda x: " -I ".join([i for i in x[0] if i.endswith("vcf")])),
          "all_stats_input": self.M2.get_output("scatter_vcf", lambda x: " -stats ".join([i for i in x[0] if i.endswith("vcf.stats")]))
        },
        script=[
          "set -euxo pipefail",
          """
        gatkm="/gatk/gatk --java-options -Xmx4g"
        echo "-----------------------learn read orientation model-----------------------"
        $gatkm LearnReadOrientationModel \
            -I ${all_f1r2_input} \
            -O artifact-priors.tar.gz 
        echo "----------------------- merge vcfs----------------------------"
        $gatkm MergeVcfs \
            -I ${all_vcf_input} \
            -O merged_unfiltered.vcf
        echo "------------------------merge mutect stats-------------------------------"
        $gatkm MergeMutectStats \
            -stats ${all_stats_input} \
            -O merged.stats
        echo --------------filter mutect stats ----------------
        $gatkm FilterMutectCalls \
            -V merged_unfiltered.vcf \
            -R ${ref_fasta} \
            --contamination-table ${contamination_table} \
            --tumor-segmentation ${segments_table} \
            --ob-priors artifact-priors.tar.gz \
            -stats merged.stats \
            --filtering-stats filtering.stats \
            -O merged_filtered.vcf 
        
        """
        ],
        outputs={
          "merged_filtered_vcf": "merged_filtered.vcf*"
        },
        docker = GATK_docker_image,
        dependencies=[self.M2, self.calc_contam]
        )
      
      self.realignment = Task(
        name = "realignment",
        inputs = {
          "command_mem": 2, 
          "bam": t_bam_cloud_path,
          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
          "input_vcf": self.merge_M2.get_output("merged_filtered_vcf", lambda x: [i for i in x[0] if i.endswith("vcf")][0]),
          "realignment_index_bundle": ref_path + config["Mutect2.realignment_index_bundle"],
          "user_project" : user_project
        },
        overrides={
          "bam" : None,
          "input_vcf": None
        },
        script=[
          """
          set -euxo pipefail
          export CLOUDSDK_CONFIG=/etc/gcloud
          /gatk/gatk --java-options "-Xmx${command_mem}g" FilterAlignmentArtifacts \
            -V ${input_vcf} -R ${ref_fasta} \
            -I ${bam} --gcs-project-for-requester-pays ${user_project} \
            --bwa-mem-index-image ${realignment_index_bundle} \
            -O realigned_filtered.vcf
          """
        ],
        outputs={
          "realigned_vcf": "realigned_filtered.vcf*"
        },
        dependencies=self.merge_M2,
        docker=GATK_docker_image
      )
      
      
      
      self.funcotate = Task(
        name = "funcotate",
        inputs={
          "merged_filtered_vcf": self.merge_M2.get_output("merged_filtered_vcf", lambda x: [i for i in x[0] if i.endswith("vcf")]),
          "data_source_folder":ref_path + config["onco_folder"],
          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
          "interval_list": ref_path + config["Mutect2.intervals"],
          "transcript_selection": ref_path + config["Mutect2.funco_transcript_selection_list"]
        },
        script=[
          """
          gatkm="/gatk/gatk --java-options -Xmx10g"
          echo "-----------------------annotate-----------------------"
          $gatkm Funcotator \
            --data-sources-path ${data_source_folder} \
            --ref-version hg19 \
            --output-file-format MAF \
            -R ${ref_fasta} \
            -V ${merged_filtered_vcf} \
            -O annot_merged_filtered.maf \
            -L ${interval_list} \
             --transcript-list ${transcript_selection} \
            --remove-filtered-variants true
          """
        ],
        overrides={
          "merged_filtered_vcf":None
        },
        outputs={
          "annot":"annot_merged_filtered.maf"
        },
        docker = GATK_docker_image,
        dependencies=self.realignment
      )
        
        
