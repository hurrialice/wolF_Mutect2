
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
                     **kwargs):

    ref_files = m2_ref_files[ref_build]
    ref_files = {**ref_files, **ref_files[sequencing_type]}
    ref_files = {**ref_files, **kwargs}

    results = {}

    # Split intervals
    intervals = SplitIntervals(inputs=dict(ref_fasta=ref_files["fasta"],
                                           interval_list=ref_files["split_intervals"], 
                                           scatter_count=scatter_count))["subintervals"]

    # MuTect2 scatter          
    m2_outputs = Mutect2(inputs=dict(pair_name=pair_name,
                                     case_name=t_name,
                                     ctrl_name=n_name,
                                     t_bam=t_bam, 
                                     t_bai=t_bai,
                                     n_bam=n_bam, 
                                     n_bai=n_bai,
                                     ref_fasta=ref_files["fasta"],
                                     gnomad_vcf=ref_files["gnomad_vcf"],
                                     pon_vcf=ref_files["pon_vcf"],
                                     interval=intervals)
                         )

    # MuTect2 gather
    m2g_outputs = GatherMutect2(inputs=dict(all_vcf_input=[m2_outputs["scatter_vcf"]]
                                           )
                               )
    results["merged_unfiltered_vcf"] = m2g_outputs["merged_unfiltered_vcf"]

    # Pileup summary scatter

    # Pileup summary gather

    # Compute contamination

    # Filter variant calls

    # Filter alignment artifacts

    # Funcotator 

    return results


