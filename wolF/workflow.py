
from .reference_files import m2_ref_files
from .tasks import *

"""
    Run a MuTect2 somatic variant calling workflow on 
    a tumor-normal pair:
        * Scatter/gather MuTect2 across intervals of the genome
        * Compute contamination
        * Filter the calls: PoN filter, gnomAD, alignment artifacts, etc.
        * Annotate the calls with Funcotator

    kwargs:
        ref_build {"hg38", "hg19"}: default "hg38"
        sequencing_type {"WGS", "WES"}: default "WGS"

    Default reference files can be overridden by specifying
    them as kwargs. Important ones:
        * ref_fasta
        * gnomad_vcf
"""
def mutect2_workflow(pair_name, t_name, n_name,
                     t_bam, t_bai,
                     n_bam, n_bai,
                     ref_build="hg38",
                     sequencing_type="WGS",
                     scatter_count=10,
                     ref_files_override={},
                     ):

    # Consolidate reference files, prioritizing the
    # "override" ones
    ref_files = m2_ref_files[ref_build]
    ref_files = {**ref_files, **ref_files[sequencing_type]}
    ref_files = {**ref_files, **ref_files_override}

    results = {}

    # Split intervals
    intervals = SplitIntervals(name="M2_SplitIntervals",
                               inputs=dict(ref_fasta=ref_files["fasta"],
                                           ref_fasta_index=ref_files["fasta_idx"],
                                           ref_fasta_dict=ref_files["fasta_dict"],
                                           interval_list=ref_files["split_intervals"], 
                                           scatter_count=scatter_count
                                          )
                              )["subintervals"]
    split_labels = [str(i) for i in range(scatter_count)]

    # Get the sample names from the BAMs
    sample_names = GetSampleNames(name="M2_GetSampleNames",
                                  inputs=dict(tumor_bam=t_bam,
                                              normal_bam=n_bam
                                             )
                                 )

    # MuTect2 scatter          
    m2_outputs = Mutect2(name="M2_Mutect2",
                         inputs=dict(case_name=sample_names["tumor_name"],
                                     ctrl_name=sample_names["normal_name"],
                                     t_bam=t_bam, 
                                     t_bai=t_bai,
                                     n_bam=n_bam, 
                                     n_bai=n_bai,
                                     ref_fasta=ref_files["fasta"],
                                     ref_fasta_index=ref_files["fasta_idx"],
                                     ref_fasta_dict=ref_files["fasta_dict"],
                                     gnomad_vcf=ref_files["gnomad_vcf"],
                                     gnomad_vcf_idx=ref_files["gnomad_vcf_idx"],
                                     pon_vcf=ref_files["pon_vcf"],
                                     pon_vcf_idx=ref_files["pon_vcf_idx"],
                                     interval=intervals,
                                     split_label=split_labels
                                    )
                         )

    # gather MuTect2 VCFs
    m2g_vcf = MergeVCFs(name="M2_MergeVCFs",
                        inputs={
                  "all_vcfs": [m2_outputs["scatter_vcf"]],
                  }
              )["merged_unfiltered_vcf"]
    results["unfiltered_vcf"] = m2g_vcf

    # gather MuTect2 callstats
    m2g_stats = MergeMutectStats(name="M2_MergeMutectStats",
                                 inputs={
                    "all_stats": [m2_outputs["scatter_vcf_stats"]],
                    }
                )["merged_stats"]
    results["unfiltered_callstats"] = m2g_stats

    # Unpack reference build
    ref_fasta = ref_files["fasta"]
    ref_fasta_idx = ref_files["fasta_idx"]
    ref_fasta_dict = ref_files["fasta_dict"]

    # Pileup summary scatter (for tumor and normal, separately)
    tumor_pileups = GetPileupSummaries(name="M2_TumorPileupScatter",
                                       inputs={
                        "bam": t_bam, "bai": t_bai,
                        "ref_fasta": ref_fasta,
                        "ref_fasta_index": ref_fasta_idx,
                        "ref_fasta_dict": ref_fasta_dict,
                        "contamination_vcf": ref_files["contamination_vcf"],
                        "contamination_vcf_idx": ref_files["contamination_vcf_idx"],
                        "interval" : intervals,
                        "split_label": split_labels,
                        "command_mem": "4",
                        }
                    )["pileups"]
    normal_pileups = GetPileupSummaries(name="M2_NormalPileupScatter",
                                        inputs={
                        "bam": n_bam, "bai": n_bai,
                        "ref_fasta": ref_fasta,
                        "ref_fasta_index": ref_fasta_idx,
                        "ref_fasta_dict": ref_fasta_dict,
                        "contamination_vcf": ref_files["contamination_vcf"],
                        "contamination_vcf_idx": ref_files["contamination_vcf_idx"],
                        "interval" : intervals,
                        "split_label": split_labels,
                        "command_mem": "4",
                        }
                    )["pileups"]

    # Pileup summary gather
    tumor_pileup = GatherPileupSummaries(name="M2_TumorPileupGather",
                                         inputs={
                       "all_pileups": [tumor_pileups],
                       "ref_fasta_dict": ref_fasta_dict,
                       }
                   )["gathered_pileup"]
    normal_pileup = GatherPileupSummaries(name="M2_NormalPileupGather",
                                          inputs={
                       "all_pileups": [normal_pileups],
                       "ref_fasta_dict": ref_fasta_dict,
                       }
                   )["gathered_pileup"]

    # Calculate contamination
    contam = CalculateContamination(name="M2_CalculateContamination",
                                    inputs={
                 "tumor_pile_tsv": tumor_pileup,
                 "normal_pile_tsv": normal_pileup
                 }
             )

    # Learn the read orientation model:
    ro_model = GatherLearnReadOrientationModel(name="M2_LearnReadOrientationModel",
                                               inputs={
                   "all_f1r2_input": [m2_outputs["f1r2"]]
                   }
               )

    # Filter variant calls
    filtered = FilterMutect2(name="M2_FilterMutect2",
                             inputs={
                   "vcf": m2g_vcf,
                   "vcf_stats": m2g_stats,
                   "ref_fasta": ref_fasta,
                   "ref_fasta_idx": ref_fasta_idx,
                   "ref_fasta_dict": ref_fasta_dict,
                   "contamination_table": contam["contamination_table"],
                   "segments_table": contam["segments_table"],
                   "artifact_priors_targz": ro_model["artifact_priors_targz"] 
                   }
               )
    results["filter_stats"] = filtered["filter_stats"]
    results["filtered_vcf"] = filtered["filtered_vcf"]

    # Filter alignment artifacts
    # LEAVE THIS OUT FOR NOW --
    # TODO REVISIT LATER
    #filtered = FilterAlignmentArtifacts()

    # Funcotator 
    func = Funcotator(name="M2_Funcotator",
                      inputs={
               "vcf": filtered["filtered_vcf"],
               "data_sources_dir": ref_files["funco_data_sources_dir"],
               "transcript_selection": ref_files["funco_transcript_selection_list"],
               "ref_build": ref_build,
               "ref_fasta": ref_fasta,
               "ref_fasta_idx": ref_fasta_idx,
               "ref_fasta_dict": ref_fasta_dict,
               "output_format": "MAF",
               }
           )
    results["funcotator_maf"] = func["annotated_output"]

    return results


