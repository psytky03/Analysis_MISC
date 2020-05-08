# README

Last modified: 2019-11-28

## Purpose

Re-align unmapped WGS reads to two HHV-6 reference genomes (Z29 and U1102) and report coverage, depth for whole viral genome and DR, U regions.  


## Prerequisite 
Please first supply the path for samtools, bwa, bcftools and name of output folder in align2HHV6_v3.sh


## Usage: run test data

`source align2HHV6_v3.sh`

Example 1: run alignment to HHV6 genome with unmapped bam from WGS

`align2hhv6_unmap_bam TestData/test1_unmap.bam test1 > test1_report.txt`


output in test1_report.txt 

```text
subtype=U1102 cov=0.93517 dp=16.6099 cov_dr=0.663865 dp_dr=14.2137 cov_u=0.899675 dp_u=15.8499 n_var=7573 file=HHV6/test1.U1102.bam
subtype=Z29 cov=0.992631 dp=18.0456 cov_dr=0.952235 dp_dr=19.4749 cov_u=0.937207 dp_u=16.9226 n_var=937 file=HHV6/test1.Z29.bam
```


Example 2: run alignment to HHV6 genome with unmapped pair-end fastq files

`align2hhv6_unmap_fq_pe TestData/test2_unmap_1.fastq TestData/test2_unmap_2.fastq test2 > test2_report.txt`

output in test1_report.txt 
subtype=U1102 cov=0.741506 dp=1.5476 cov_dr=0.512301 dp_dr=1.2763 cov_u=0.714115 dp_u=1.47936 n_var=5091 file=HHV6/test2.U1102.bam
subtype=Z29 cov=0.796835 dp=1.69786 cov_dr=0.768793 dp_dr=2.05448 cov_u=0.752632 dp_u=1.57974 n_var=658 file=HHV6/test2.Z29.bam

