
import json
import pandas as pd

from wolf import Task, ReadFile

GATK_docker_image = 'gcr.io/broad-getzlab-workflows/gatk4_wolf:v6'

class SplitIntervals(Task):
    name = "SplitIntervals"

    inputs = {
      "ref_fasta": None,
      "ref_fasta_index": None,
      "ref_fasta_dict": None,
      "interval_list" : None, # .list, .intervals, .bed, or .vcf
      "scatter_count" : 10,
    }

    script = """
        set -euxo pipefail

        ln -s ${ref_fasta} refbuild.fa
        ln -s ${ref_fasta_index} refbuild.fa.fai
        ln -s ${ref_fasta_dict} refbuild.dict

        gatk SplitIntervals -R refbuild.fa -L ${interval_list} \
          --scatter-count ${scatter_count} -O intervals
         """
    
    outputs = {
      "subintervals" : "intervals/*-scattered.interval_list"
    }

    docker = GATK_docker_image


class GetSampleNames(Task):
    name="GATK_GetSampleNames"

    inputs={"tumor_bam" : None,
            "normal_bam" : None
    }

    script = """
        set -euxo pipefail
        gatk GetSampleName -I ${normal_bam} -O normal_name.txt 
        gatk GetSampleName -I ${tumor_bam} -O tumor_name.txt 
    """

    outputs = {"normal_name": ("normal_name.txt", ReadFile),
               "tumor_name" : ("tumor_name.txt", ReadFile)
               }
    
    docker = GATK_docker_image


class Mutect2(Task):
    name="Mutect2"
    inputs={
        "case_name": None,
        "ctrl_name": None,
        "t_bam": None, "t_bai": None,
        "n_bam": None, "n_bai": None,
        "ref_fasta": None,
        "ref_fasta_index": None,
        "ref_fasta_dict": None,
        "gnomad_vcf": None,
        "gnomad_vcf_idx": None,
        "pon_vcf": None,
        "pon_vcf_idx": None,
        "command_mem": "7",
        "interval" : "",
        "split_label": "0",
        "extra_args": "",
    }

    script = """ 
        set -euxo pipefail

        # Figure out whether a set of intervals was provided
        if test ${interval}
        then
            interval_str="-L ${interval}"
        else
            interval_str=""
        fi
        
        # Make sure reference fasta/fai/dict are set up correctly
        ln -s ${ref_fasta} refbuild.fa
        ln -s ${ref_fasta_index} refbuild.fa.fai
        ln -s ${ref_fasta_dict} refbuild.dict

        # Make sure the bams/bais are named correctly
        ln -s ${t_bam} tumor.bam
        ln -s ${t_bai} tumor.bam.bai
        ln -s ${n_bam} normal.bam
        ln -s ${n_bai} normal.bam.bai

        # Make sure PoN VCF/index are set up correctly
        pon_ext=${pon_vcf#*.}
        pon_path="pon.${pon_ext}"
        ln -s ${pon_vcf} ${pon_path}
        
        pon_idx_ext=${pon_vcf_idx#*.}
        ln -s ${pon_vcf_idx} pon.${pon_idx_ext}

        # Make sure the gnomAD VCF/index are set up correctly
        gnomad_ext=${gnomad_vcf#*.}
        gnomad_path="gnomad.${gnomad_ext}"
        ln -s ${gnomad_vcf} ${gnomad_path}

        gnomad_idx_ext=${gnomad_vcf_idx#*.}
        ln -s ${gnomad_vcf_idx} gnomad.${gnomad_idx_ext}

        echo ${interval}
        gatk --java-options "-Xmx${command_mem}g" Mutect2 \
            -R refbuild.fa \
            -I tumor.bam -tumor ${case_name} \
            -I normal.bam -normal ${ctrl_name} \
            --germline-resource ${gnomad_path} \
            -pon ${pon_path} \
            ${interval_str} \
            -O scatter_${split_label}.vcf \
            --f1r2-tar-gz f1r2_${split_label}.tar.gz \
            ${extra_args}
    """

    outputs = {
      "scatter_vcf": "scatter_*.vcf",
      "scatter_vcf_stats": "*.stats",
      "f1r2": "f1r2_*.tar.gz"
    }
    
    resources = {
        "mem" : "8G",
        "cpus-per-task" : 1
    }
    docker = GATK_docker_image


class MergeVCFs(Task):
   
    inputs = {"all_vcfs": None}

    script = """
        set -euxo
        
        # Need all vcf files, separated by " -I "
        vcfs_str=$(python -c 'import sys;print(" -I ".join([l.strip() for l in sys.stdin]))' < ${all_vcfs})

        gatk MergeVcfs \
            -I $vcfs_str \
            -O merged_unfiltered.vcf
    """

    outputs = {"merged_unfiltered_vcf": "merged_unfiltered.vcf"
              }

    resources = {"mem": "4GB"
                }

    docker = GATK_docker_image


class MergeMutectStats(Task):
   
    inputs = {"all_stats": None}

    script = """
        set -euxo

        # Need all stats files, separated by " -stats "
        stats_str=$(python -c 'import sys;print(" -stats ".join([l.strip() for l in sys.stdin]))' < ${all_stats})

        gatk MergeMutectStats \
            -stats $stats_str \
            -O merged.stats
    """

    outputs = {"merged_stats": "merged.stats"}

    docker = GATK_docker_image



class GetPileupSummaries(Task):
    
    inputs = {
        "bam": None, "bai": None,
        "ref_fasta": None,
        "ref_fasta_index": None,
        "ref_fasta_dict": None,
        "contamination_vcf": None,
        "contamination_vcf_idx": None,
        "interval" : "",
        "command_mem": "4",
        "split_label": "0"
    }

    script = """
        set -euxo
        
        # Figure out whether a set of intervals was provided
        if test ${interval}
        then
            interval_str="-L ${interval}"
        else
            interval_str=""
        fi

        # Enforce BAM/BAI naming conventions
        ln -s ${bam} input.bam
        ln -s ${bai} input.bai

        # Enforce FASTA naming conventions
        ln -s ${ref_fasta} ref.fa
        ln -s ${ref_fasta_index} ref.fa.fai
        ln -s ${ref_fasta_dict} ref.dict

        # Enforce VCF index naming conventions
        contam_ext=${contamination_vcf#*.}
        contamination_vcf_path="contam.$contam_ext"
        ln -s ${contamination_vcf} ${contamination_vcf_path}

        contam_idx_ext=${contamination_vcf_idx#*.}
        ln -s ${contamination_vcf_idx} contam.${contam_idx_ext}

        # Finally: run GATK GetPileupSummaries
        gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries -R ref.fa \
          -I ${bam} --interval-set-rule INTERSECTION \
          ${interval_str} \
          -V ${contamination_vcf_path} -L ${contamination_vcf_path} \
          -O pileups_${split_label}.table
    """

    outputs = {"pileups": "pileups_*.table"
              }

    resources = {"mem": "4GB"
                }

    docker = GATK_docker_image


class GatherPileupSummaries(Task):

    inputs = {"all_pileups": None,
              "ref_fasta_dict": None
             }

    script = """
        set -euxo

        # Need all pileup files, separated by "-I"
        pileups_str=$(python -c 'import sys;print(" -I ".join([l.strip() for l in sys.stdin]))' < ${all_pileups})

        gatk --java-options -Xmx2g GatherPileupSummaries \
                -I $pileups_str \
                --sequence-dictionary ${ref_fasta_dict} \
                -O pile.tsv
    """


    outputs = {"gathered_pileup": "pile.tsv"
              }

    docker = GATK_docker_image


class CalculateContamination(Task):

    inputs = {"tumor_pile_tsv": None,
              "normal_pile_tsv": None,
              }

    script = """
        set -euxo
        gatk --java-options -Xmx2g CalculateContamination \
            -I ${tumor_pile_tsv} -O contamination.table \
            --tumor-segmentation segments.table \
            -matched ${normal_pile_tsv}
    """
    
    outputs = {"contamination_table": "contamination.table",
               "segments_table": "segments.table"
              }

    docker = GATK_docker_image


class GatherLearnReadOrientationModel(Task):
    
    inputs = {"all_f1r2_input": None,
              "command_mem": 15, 
             }

    script = """
        set -euxo

        # Need all f1r2 files, separated by " -I "
        f1r2_str=$(python -c 'import sys;print(" -I ".join([l.strip() for l in sys.stdin]))' < ${all_f1r2_input})

        gatk LearnReadOrientationModel \
            --java-options "-Xmx${command_mem}g" \
            -I $f1r2_str \
            -O artifact-priors.tar.gz 
    """

    outputs = {"artifact_priors_targz": "artifact-priors.tar.gz"
              }

    resources = { "mem": "16GB"
                }

    docker = GATK_docker_image


class FilterMutect2(Task):

    inputs = {"vcf": None,
              "vcf_stats": None,
              "ref_fasta": None,
              "ref_fasta_idx": None,
              "ref_fasta_dict": None,
              "contamination_table": None,
              "segments_table": None,
              "artifact_priors_targz": None,
              "command_mem": 7
             }

    script = """
        set -euxo

        # Standardize ref fasta file names
        ln -s ${ref_fasta} refbuild.fa
        ln -s ${ref_fasta_idx} refbuild.fa.fai
        ln -s ${ref_fasta_dict} refbuild.dict

        gatk FilterMutectCalls \
            --java-options "-Xmx${command_mem}g" \
            -V ${vcf} \
            -R refbuild.fa \
            --contamination-table ${contamination_table} \
            --tumor-segmentation ${segments_table} \
            --ob-priors ${artifact_priors_targz} \
            -stats ${vcf_stats} \
            --filtering-stats filtering.stats \
            -O merged_filtered.vcf 
    """
    
    outputs = {"filtered_vcf": "merged_filtered.vcf",
               "filter_stats": "filtering.stats"
              }
   
    resources = {"mem": "8GB"
                }

    docker = GATK_docker_image


class FilterAlignmentArtifacts(Task):
    #TODO FINISH IMPLEMENTING THIS 
    inputs = {"vcf": None,
              "tumor_bam": None,
              "tumor_bai": None,
              "ref_fasta": None,
              "ref_fasta_idx": None,
              "ref_fasta_dict": None,
              "bwamem_index_image": None,
              "command_mem": 7
             }

    script = """
        set -euxo pipefail

        gatk --java-options "-Xmx${command_mem}g" FilterAlignmentArtifacts \
          -V ${vcf} -R refbuild.fa \
          -I ${bam} \
          --bwa-mem-index-image ${bwamem_index_image} \
          -O realigned_filtered.vcf

    """
    
    outputs = {"filtered_vcf": "realigned_filtered.vcf"}

    resources = {"mem": "8GB"
                }

    docker = GATK_docker_image


class Funcotator(Task):

    inputs = {"vcf": None,
              "data_sources_dir": None,
              "transcript_selection": None,
              "ref_build": None,
              "ref_fasta": None,
              "ref_fasta_idx": None,
              "ref_fasta_dict": None,
              "command_mem": 15,
              "output_format": "MAF",
            }

    script = """
        set -euxo

        # Standardize fasta filenames
        ln -s ${ref_fasta} ref.fa
        ln -s ${ref_fasta_idx} ref.fa.fai
        ln -s ${ref_fasta_dict} ref.dict

        gatk --java-options "-Xmx${command_mem}g" Funcotator \
          --data-sources-path ${data_sources_dir} \
          --transcript-list ${transcript_selection} \
          --ref-version ${ref_build} \
          --output-file-format ${output_format} \
          --remove-filtered-variants true \
          -V ${vcf} \
          -R ref.fa \
          -O annot.${output_format} \
    """
    
    outputs = {"annotated_output": "annot.*"
              }

    resources = {"mem": "16GB"
                }

    docker = GATK_docker_image


############################################################
# OLD
#############################################################
#
## read the input.json file for configuration
#class Mutect2(Workflow):
#
#    def workflow(self, 
#      t_bam_cloud_path, n_bam_cloud_path,
#      user_project):
#      
#      with open('mutect2.inputs.json', 'r') as f:
#        config  = json.load(f)
#      
#      ref_path = config["ref_local_path"]
#      
#      
#      self.split_intervals = Task(
#        name = "split_intervals",
#        inputs = {
#          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
#          "scatter_count" : config['Mutect2.scatter_count'],
#          "wes_intervals" : ref_path + config['Mutect2.intervals'],
#          "user_project": user_project
#        },
#        script = [
#          "set -euxo pipefail",
#          "export CLOUDSDK_CONFIG=/etc/gcloud",
#          """/gatk/gatk SplitIntervals -R ${ref_fasta} -L ${wes_intervals} \
#            --scatter-count ${scatter_count} -O intervals"""
#        ],
#        outputs = {
#          "subintervals" : "intervals/*-scattered.interval_list"
#        },
#        overrides={
#          "wes_intervals": None,
#          "ref_fasta" : None
#        },
#        docker = GATK_docker_image
#      )
#      
#      self.get_sample_name = Task(
#        name = "get_sample_name",
#        inputs = {
#          "tumor_bam" : t_bam_cloud_path,
#          "normal_bam" : n_bam_cloud_path,
#          "user_project" : user_project
#        },
#        script = [
#          "set -euxo pipefail",
#          "export CLOUDSDK_CONFIG=/etc/gcloud",
#          "/gatk/gatk GetSampleName -I ${normal_bam} -O normal_name.txt --gcs-project-for-requester-pays ${user_project}",
#          "/gatk/gatk GetSampleName -I ${tumor_bam} -O tumor_name.txt --gcs-project-for-requester-pays ${user_project} "
#        ],
#        outputs = {
#          "tumor_name" : ("tumor_name.txt", output_helpers.read_file),
#          "normal_name" : ("normal_name.txt", output_helpers.read_file)
#        },
#        overrides={
#          "tumor_bam": None,
#          "normal_bam": None,
#          "user_project": None
#        },
#        docker = GATK_docker_image
#      )
#      
#      
#      self.M2 = Task(
#        name="scatter_m2",
#        inputs={
#          "ref_fasta": ref_path + config['Mutect2.ref_fasta'],
#          "gnomad": ref_path + config['Mutect2.gnomad'],
#          "default_pon": ref_path + config['Mutect2.pon'],
#          "variants_for_contamination": ref_path + config['Mutect2.variants_for_contamination'],
#          "tumor_bam": t_bam_cloud_path,
#          "normal_bam": n_bam_cloud_path,
#          "tumor_name": self.get_sample_name.get_output("tumor_name"),
#          "normal_name": self.get_sample_name.get_output("normal_name"),
#          "extra_args": config['Mutect2.m2_extra_args'],
#          "user_project": user_project,
#          "command_mem": "4",
#          "interval" : self.split_intervals.get_output("subintervals", lambda x: x[0])
#        },
#        script = [
#          "set -euxo pipefail",
#          "export CLOUDSDK_CONFIG=/etc/gcloud",
#          "echo $interval",
#          """
#          /gatk/gatk --java-options "-Xmx${command_mem}g" Mutect2 \
#            -R ${ref_fasta} \
#            -I ${tumor_bam} -tumor ${tumor_name} \
#            -I ${normal_bam} -normal ${normal_name} \
#            --germline-resource ${gnomad} \
#            -pon ${default_pon} \
#            -L  ${interval} \
#            -O scatter.vcf \
#            --f1r2-tar-gz f1r2.tar.gz \
#            --gcs-project-for-requester-pays ${user_project}\
#            ${extra_args}
#          """,
#          
#          """
#          /gatk/gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries -R ${ref_fasta}\
#            -I ${tumor_bam} --interval-set-rule INTERSECTION -L ${interval} \
#            -V ${variants_for_contamination} -L ${variants_for_contamination} \
#            --gcs-project-for-requester-pays ${user_project}\
#            -O tumor-pileups.table
#          """,
#          
#          """
#          /gatk/gatk --java-options "-Xmx${command_mem}g" GetPileupSummaries -R ${ref_fasta}\
#            -I ${normal_bam} --interval-set-rule INTERSECTION -L ${interval} \
#            -V ${variants_for_contamination} -L ${variants_for_contamination}\
#            --gcs-project-for-requester-pays ${user_project}\
#            -O normal-pileups.table
#          """
#        ],
#        outputs = {
#          "scatter_vcf": "scatter.vcf*",
#          "tpile": "tumor-pileups.table",
#          "npile": "normal-pileups.table",
#          "f1r2": "f1r2.tar.gz"
#        },
#        
#        resources = {
#            "mem" : "6G",
#            "cpus-per-task" : 1
#        },
#        overrides={
#          "normal_bam": None,
#          "tumor_bam": None
#        },
#        docker = GATK_docker_image,
#        dependencies = [self.split_intervals, self.get_sample_name]
#      )
#      self.calc_contam = Task(
#        name = "calculate_contamination",
#        inputs={
#          "all_tumor_pile_input": self.M2.get_output("tpile", lambda x: " -I ".join(x)),
#          "all_normal_pile_input": self.M2.get_output("npile", lambda x: " -I ".join(x)),
#          "seq_dict": ref_path + config['Mutect2.ref_dict']
#        },
#        script=[
#          "set -euxo pipefail",
#          """
#          /gatk/gatk --java-options -Xmx2g GatherPileupSummaries -I ${all_tumor_pile_input} \
#            --sequence-dictionary ${seq_dict} -O tumor_pile.tsv
#          """,
#          """
#          /gatk/gatk --java-options -Xmx2g GatherPileupSummaries -I ${all_normal_pile_input} \
#            --sequence-dictionary ${seq_dict} -O normal_pile.tsv
#          """,
#          """
#          /gatk/gatk --java-options -Xmx2g CalculateContamination \
#            -I tumor_pile.tsv -O contamination.table \
#            --tumor-segmentation segments.table \
#            -matched normal_pile.tsv
#          """
#        ],
#        outputs={
#          "contamination": "contamination.table",
#          "segments": "segments.table",
#          "tumor_pile_table" : "tumor_pile.tsv",
#          "normal_pile_table" : "normal_pile.tsv"
#        },
#        docker = GATK_docker_image,
#        dependencies=self.M2)
#      
#      self.merge_M2 = Task(
#        name="merge_m2",
#        inputs={
#          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
#          "contamination_table": self.calc_contam.get_output("contamination"),
#          "segments_table": self.calc_contam.get_output("segments"),
#          "all_f1r2_input": self.M2.get_output("f1r2", lambda x: " -I ".join(x)),
#          "all_vcf_input": self.M2.get_output("scatter_vcf", lambda x: " -I ".join([i for i in x[0] if i.endswith("vcf")])),
#          "all_stats_input": self.M2.get_output("scatter_vcf", lambda x: " -stats ".join([i for i in x[0] if i.endswith("vcf.stats")]))
#        },
#        script=[
#          "set -euxo pipefail",
#          """
#        gatkm="/gatk/gatk --java-options -Xmx4g"
#        echo "-----------------------learn read orientation model-----------------------"
#        $gatkm LearnReadOrientationModel \
#            -I ${all_f1r2_input} \
#            -O artifact-priors.tar.gz 
#        echo "----------------------- merge vcfs----------------------------"
#        $gatkm MergeVcfs \
#            -I ${all_vcf_input} \
#            -O merged_unfiltered.vcf
#        echo "------------------------merge mutect stats-------------------------------"
#        $gatkm MergeMutectStats \
#            -stats ${all_stats_input} \
#            -O merged.stats
#        echo --------------filter mutect stats ----------------
#        $gatkm FilterMutectCalls \
#            -V merged_unfiltered.vcf \
#            -R ${ref_fasta} \
#            --contamination-table ${contamination_table} \
#            --tumor-segmentation ${segments_table} \
#            --ob-priors artifact-priors.tar.gz \
#            -stats merged.stats \
#            --filtering-stats filtering.stats \
#            -O merged_filtered.vcf 
#        
#        """
#        ],
#        outputs={
#          "merged_filtered_vcf": "merged_filtered.vcf*"
#        },
#        docker = GATK_docker_image,
#        dependencies=[self.M2, self.calc_contam]
#        )
#      
#      self.realignment = Task(
#        name = "realignment",
#        inputs = {
#          "command_mem": 2, 
#          "bam": t_bam_cloud_path,
#          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
#          "input_vcf": self.merge_M2.get_output("merged_filtered_vcf", lambda x: [i for i in x[0] if i.endswith("vcf")][0]),
#          "realignment_index_bundle": ref_path + config["Mutect2.realignment_index_bundle"],
#          "user_project" : user_project
#        },
#        overrides={
#          "bam" : None,
#          "input_vcf": None
#        },
#        script=[
#          """
#          set -euxo pipefail
#          export CLOUDSDK_CONFIG=/etc/gcloud
#          /gatk/gatk --java-options "-Xmx${command_mem}g" FilterAlignmentArtifacts \
#            -V ${input_vcf} -R ${ref_fasta} \
#            -I ${bam} --gcs-project-for-requester-pays ${user_project} \
#            --bwa-mem-index-image ${realignment_index_bundle} \
#            -O realigned_filtered.vcf
#          """
#        ],
#        outputs={
#          "realigned_vcf": "realigned_filtered.vcf*"
#        },
#        dependencies=self.merge_M2,
#        docker=GATK_docker_image
#      )
#      
#      
#      
#      self.funcotate = Task(
#        name = "funcotate",
#        inputs={
#          "merged_filtered_vcf": self.merge_M2.get_output("merged_filtered_vcf", lambda x: [i for i in x[0] if i.endswith("vcf")]),
#          "data_source_folder":ref_path + config["onco_folder"],
#          "ref_fasta": ref_path + config["Mutect2.ref_fasta"],
#          "interval_list": ref_path + config["Mutect2.intervals"],
#          "transcript_selection": ref_path + config["Mutect2.funco_transcript_selection_list"]
#        },
#        script=[
#          """
#          gatkm="/gatk/gatk --java-options -Xmx10g"
#          echo "-----------------------annotate-----------------------"
#          $gatkm Funcotator \
#            --data-sources-path ${data_source_folder} \
#            --ref-version hg19 \
#            --output-file-format MAF \
#            -R ${ref_fasta} \
#            -V ${merged_filtered_vcf} \
#            -O annot_merged_filtered.maf \
#            -L ${interval_list} \
#             --transcript-list ${transcript_selection} \
#            --remove-filtered-variants true
#          """
#        ],
#        overrides={
#          "merged_filtered_vcf":None
#        },
#        outputs={
#          "annot":"annot_merged_filtered.maf"
#        },
#        docker = GATK_docker_image,
#        dependencies=self.realignment
#      )
#        
#        
