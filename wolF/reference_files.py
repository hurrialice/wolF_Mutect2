
m2_ref_files = {
    "hg38": {
        "fasta"     : "gs://getzlab-workflows-reference_files-oa/hg38/gdc/GRCh38.d1.vd1.fa",
        "fasta_idx" : "gs://getzlab-workflows-reference_files-oa/hg38/gdc/GRCh38.d1.vd1.fa.fai",
        "fasta_dict": "gs://getzlab-workflows-reference_files-oa/hg38/gdc/GRCh38.d1.vd1.dict",
        
        # TODO: gnomAD seems old, ca. 2017
        "gnomad_vcf": "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz",
        "gnomad_vcf_idx": "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi",
        "pon_vcf": "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz",
        "pon_vcf": "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi",
        "contamination_vcf": "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz",
        "contamination_vcf_idx": "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi",

        # Funcotator
        "funco_data_sources_dir"    : "gs://getzlab-workflows-reference_files-oa/funcotator/funcotator_dataSources.v1.7.20200521s/",
        "funco_transcript_selection_list" :  "gs://broad-public-datasets/funcotator/transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt",

        "WGS": {
            "split_intervals": "gs://getzlab-workflows-reference_files-oa/hg38/final_hg38_wgs_masked_intervals.bed",
        },
        "WES": {
            "split_intervals": "gs://gcp-public-data--broad-references/hg38/v0/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list",
        }

    },
    "hg19": {
        "ref_fasta"    : "gs://getzlab-workflows-reference_files-oa/hg19/Homo_sapiens_assembly19.fasta",
        "ref_fasta_idx" : "gs://getzlab-workflows-reference_files-oa/hg19/Homo_sapiens_assembly19.fasta.fai",
        "ref_fasta_dict": "gs://getzlab-workflows-reference_files-oa/hg19/Homo_sapiens_assembly19.dict",
       
        # TODO: gnomAD seems old, ca. 2017
        "gnomad_vcf": "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf",
        "gnomad_vcf_idx": "gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx",
        "contamination_vcf": "gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf",
        "contamination_vcf_idx": "gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf.idx",
        
        # Funcotator
        "funco_data_sources_dir"    : "gs://getzlab-workflows-reference_files-oa/funcotator/funcotator_dataSources.v1.7.20200521s/",
        "funco_transcript_selection_list" :  "gs://broad-public-datasets/funcotator/transcriptList.exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt",

        "WGS": {
            "split_intervals": "gs://broad-getzlab-wolf-wgs-hg19/additional_refs/M2.interval.CytoBand.acen.Blacklist.mask.12Oct2016.bed",
            "pon_vcf": "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf",
            "pon_vcf_idx": "gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx",
        },
        "WES": {
            "split_intervals": "gs://getzlab-workflows-reference_files-oa/hg19/agilent/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list",
            "pon_vcf": "gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf",
            "pon_vcf_idx": "gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx",
        }
    }
}

