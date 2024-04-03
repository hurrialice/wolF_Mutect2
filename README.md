# GATK Mutect2 in wolF

## Info

- Docker used: `gcr.io/broad-getzlab-workflows/gatk4_wolf:v6`
- Reference files - Most come from [GATK-best-practice bundle](https://console.cloud.google.com/storage/browser/gatk-best-practices/). See `wolF/reference_files.py` for details.
- Original WDL for [GATK 4.1.4.1 Best Practice](https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/tree/2.7.0)
- Genome build: hg38 *or* hg19

## News

2024-03. Reimplemented the workflow with modern wolF design patterns

## Acknowledgement

Thanks to Aaron for building robust [canine](https://github.com/broadinstitute/canine) backend, to Julian for creating [wolF](https://github.com/getzlab/wolF) and quick & responsive support.


