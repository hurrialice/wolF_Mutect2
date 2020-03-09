# GATK Mutect2 in wolF

## Info

- Docker used: broadinstitute/gatk:4.1.4.1
- Reference files - including PON, come from [GATK-best-practice bundle](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37)
- Original WDL for [GATK 4.1.4.1 Best Practice](https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/tree/2.7.0)
- Genome build: hg19
- Disk image: 
    - NAME: mutect2-wolf
    - PROJECT: broad-cga-qing-gdac 

## News

A minimal implementation tested work, will add 

- [x] Build a disk image for all reference files, ~~potentially as another task~~
- [ ] Cross-validate on Terra with same input
- [ ] Cost estimate
- [ ] Incorporates more optional arguments (bells and whistles)

## Acknowledgement

Thanks to Aaron for building robust [canine](https://github.com/broadinstitute/canine) backend, to Julian for creating [wolF](https://github.com/getzlab/wolF) and quick & responsive support.


