# Zenomix
Zenome Pharmacogenomics caller for mixed NGS data (WES + LP WGS), forked from StellarPGx


CYP450 genes supported: CYP2D6, CYP2B6, CYP2C19, CYP2C9, CYP3A4, CYP3A5, CYP1A2
Non-CYP450 genes supported: CYPOR (POR)

Zenomix is based on StellarPGx (please refer to the excellent documentation [here](https://github.com/SBIMB/StellarPGx)). Main script is built using PHP and utilize Docker technology.

The program is focused on analysis of short-read high-coverage whole exome sequence data in combination with low-coverage(2-7x) whole genome sequence data. To analyze hybrid alleles in such complicated gene as CYP2D6, we developed hybrid algorithm that extracts reads from certain parts of paralog, which has high homology to the target gene. Developed especially to analyze the metabolism of drugs used in psychiatry. 7 determined cytochromes are involved in the inactivation of more than 95% of all known psychiatric drugs. Cytochrome P450 Oxidoreductase donate electrons directly from NADPH to all microsomal P450 enzymes. Therefore, impaired POR function can alter the metabolism of all P450 cytochromes.


# Getting started
The following are required to run the Zenomix pipeline:

1. Prerequisite software
  - PHP (preferably v8.x)
  - Docker

2. Whole genome sequence (WGS) data
  - Indexed BAM/CRAM files

3. Reference genome
  - hg19, b37, or hg38


## Installation
### PHP:
Install PHP on Linux by running the following command (Skip if you have PHP installed already):
```
sudo apt-get install php
```
### Docker:
For Docker installation, please refer to the excellent documentation [here](https://docs.docker.com/get-docker/))

### Zenomix repository:
Clone the StellarPGx repository by running the following command:
```
git clone https://github.com/zenomeplatform/pharm_pipe.git && cd pharm_pipe
```

## Running Zenomix
### Step 1 - Parameters
The parameters for Docker are set as default.

### Step 2 - Execution of the pipeline
For execution on a local machine
```
php -f run_zenomix.php ref_file bam_path bam_name work_mode report_path
```
### Step 3 - Expected output
The expected output excel file contains genotypes and phenotypes for 7 CYP450 genes and CYPOR. Additional sheets contain phenotypes for 7 groups of drugs: "Antidepressants", "Anticonvulsants", "Neuroleptics", "Tranquilizers", "Sedatives", "Stimulants" and "Nootropic". Presumptive impact of POR genotype on cytochromes phenotypes is provided in "POR_Impact" sheet. Some variants (intronic, splice site, promoter) cannot be accurately identified by low-coverage(2-7x) WGS, so there is a risk of having alleles which are defined by these core variants. Information about risk is collected in "Risk_alleles" sheet.

![show_result](https://user-images.githubusercontent.com/91198710/213020784-935a5793-ec34-4f83-a6ed-01bdcb7fdabf.png)
