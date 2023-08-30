
<!-- # Badges -->
<p align="center">
    <a href="https://www.nextflow.io/">
        <img src="https://img.shields.io/badge/nextflow-23.04.2-%233EAA26" /></a>
    <a href="https://www.docker.com/">
        <img src="https://img.shields.io/badge/Docker-blue?logo=docker&logoColor=white" /></a>
    <a href="https://docs.conda.io/en/latest/">
        <img src="https://img.shields.io/badge/conda-%233eb049?logo=anaconda&logoColor=white" /></a>
    <a href="https://www.python.org/">
        <img src="https://img.shields.io/badge/python-%23306998?logo=python&logoColor=white" /></a>
</p>

<br>
<br>

<!-- # FidelityFinder logo -->
<p align="center">
  <img src="https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/ff_logo.jpg?raw=true" alt="Logo">
</p>


## Table of Contents 

1. [Introduction](#introduction)
2. [Quick Start](#quick-start)
3. [Input parameters](#input-parameters)
4. [Libray preparation](#library-preparation)
5. [Pipeline overview](#pipeline-overview)
   * [Step 1: Quality Control](#step-1-quality-control)
   * [Step 2: Joining of Paired Reads](#step-2-joining-of-paired-reads)
   * [Step 3: Graph Lengths of the Merged Reads](#step-3-graph-lengths-of-the-merged-reads)
   * [Step 4: Obtain Sequences in Fasta Format](#step-4-obtain-sequences-in-fasta-format)
   * [Step 5: Obtain Consensus Sequences](#step-5-obtain-consensus-sequences)
   * [Step 6: Map Consensus Sequences](#step-6-map-consensus-sequences)
   * [Step 7: Variant Calling](#step-7-variant-calling)
   * [Step 8: VCF Analysis](#step-8-vcf-analysis)
   * [Step 9: Offspring Search (optional)](#step-9-offspring-search-optional)
6. [Test data results](#test-data-results)
7. [Extra content](#extra-content)
8. [Creative Commons](#creative-commons)
9. [Citation](#citation)
10. [References](#references)

## Introduction
**FidelityFinder** is a bioinformatics analysis pipeline to determine the fidelity transcriptases and the fidelity of DNA synthesis by reverse transcriptases (RTs) from sequences obtained by Next Generation Sequencing (NGS). 

Fidelity determination is based on the "Primers IDs method". Each cDNA obtained by RTs is tagged with a barcode, so each cDNA molecule has a unique identity. Then, these cDNAs are amplified by PCR and added the adaptor sequences to generate a library (see [Libray preparation](#library-preparation) for a more detailed explanation). The libraries generated are sequenced then by NGS and **FidelityFinder** is able to evaluate the fidelity of the RT used: the pipeline is able to discard PCR and NGS errors thanks to the construction of consensus sequences that share the same barcode sequence, so it obtains and error rate. This error rate is the combination of transcription and reverse transcription errors (see [Pipeline overview](#pipeline-overview) for a more detailed explanation of the pipeline). The higher the fidelity of the RT used, the lower the error rate.  

Primers IDs method to determine fidelity of reverse transcriptases:

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Primers_IDs_method.PNG?raw=true)

## Quick Start 
1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. Install [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html#regular-installation)
3. Run the pipeline with the test data provided:
   ```console
   nextflow run https://github.com/friburgo-moc/FidelityFinder -profile test_profile -r main
   ```
   You can find the expected results of the test data analysis in [Test data results](#test-data-results)
   
4. Run the pipeline with your own data:
   ```console
   nextflow run https://github.com/friburgo-moc/FidelityFinder -c <config_file> -profile my_profile -r main
   ```
   `<config_file>` must contain your input parameters for your own analysis. See [Input parameters](#input-parameters) section for more details.
   
---
If you are using conda, it is recommended to use the [`conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html#how-it-works) setting to store the environments in a central location for future pipeline runs.

If you want to generate an [`Execution report`](https://www.nextflow.io/docs/latest/tracing.html#execution-report), add the following command line option: 
```console
-with-report [file name]
```

## Input parameters

Input parameters must be provided in your own config file. It must contain the following structure:

  ```console
conda.enabled = true
 my_profile  {
	params.seq_folder_path = '/folder/with/raw/sequences'
	params.ref_seq_path = 'file/with/reference/sequence'
	params.insert_length = "length of the library insert"
	params.fw_primer = "forward_primer_sequence"
	params.rv_primer = "reverse_primer_sequence"
	params.cutoff = "cutoff_value_to_discard_low_frequency_barcodes"
	params.bc_size = "nucleotides_of_the_barcode"
	params.threshold = "threshold_used_to_construct_consensus_sequences"
	params.min_pos = "first position of the reference sequence used to quantify mutations"
	params.max_pos = "last position of the reference sequence used to quantify mutations"
 }
  ```

- **params.seq_folder_path**  
path of the folder where you have the raw sequences obtained by NGS. Files must have the following names: **<sample_name>_1.fastq** and **<sample_name>_2.fastq** for the forward and reverse sequences, respectively. You can also have them compressed: **<sample_name>_1.fastq.gz** and **<sample_name>_2.fastq.gz**.
- **params.ref_seq_path**  
path of the file with the reference sequence, i.e., the sequence of the insert of the library (without mutations), the barcode sequence must be indicated with as many "N" as nucleotides it has.
- **params.insert_length**
length of the library insert that has been sequenced, i.e. total length of the library except for the adaptors. This parameter is important to filter merged reads according to their length: reads that differ by more than 20 nucleotides from the indicated length are not selected. If you want to modify the ±20 nucleotides consideration, you have to change the "20" numbers found in the following line of `main.nf`:
```console
awk 'BEGIN {FS = "\\t" ; OFS = "\\n"} {header = \$0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= ("${params.insert_length}"-20) && length(seq) <= ("${params.insert_length}"+20)) {print header, seq, qheader, qseq}}' "${sampleId}.assembled.fastq" > "${sampleId}.assembled_filtered.fastq"
```
- **params.fw_primer**  
- **params.rv_primer**  
- **params.cutoff**  
- **params.bc_size**  
- **params.threshold**
- **params.min_pos**
- **params.max_pos** 

Once you have created your config file with your own parameters, you have to run the pipeline and specify the path of the config_file with -c <config_file> and the name of the profile with your specific parameters -profile my_profile

## Library preparation

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Library_preparation_workflow.PNG?raw=true)

## Pipeline overview 

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Workflow.PNG?raw=true)

The pipeline is built using Nextflow. Processing steps:

- ### Step 1: Quality Control
  Quality analyses are performed over reads using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) software.

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `1_fastqc_results/`
    		- `<sample_name>_1_fastqc*.html`: html file
      		- `<sample_name>_1_fastqc*.zip`: zip file	 
    		- `<sample_name>_2_fastqc*.html`: html file
      		- `<sample_name>_2_fastqc*.zip`: zip file

	</details>

  
- ### Step 2: Joining of paired reads
  The first step is to join paired reads into one single sequence. This has been performed using [PEAR](https://sco.h-its.org/exelixis/web/software/pear/doc.html)
(Paired-End reAd mergeR) (J. Zhang et al., 2014), a fast and accurate paired-end read merger. PEAR evaluates all
possible paired-end read overlaps without requiring the target fragment size as input. In addition, it
implements a statistical test for minimizing false-positive results. PEAR outputs four files. A file
containing the assembled reads (assembled.fastq extension), two files containing the forward and
reverse unassembled reads (unassembled.forward.fastq and unassembled.reverse.fastq extensions),
and a file containing the discarded reads (discarded.fastq extension) (cf. Appendix 1: Checklist of
provided supplemental files, directory: pear/). We performed this method because it holds the
adapters (including the barcode), which are interesting for our analysis.

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `2_merged_sequences/`
			- `<sample_name>.unassembled.forward.fastq`:
			- `<sample_name>.unassembled.reverse.fastq`:	 
			- `<sample_name>.discarded.fast`:
			- `<sample_name>.assembled.fastq`:
			- `<sample_name>.assembled_filtered.fastq`:
			- `<sample_name>._lengths.txt`:
        
	</details>

- ### Step 3: Graph lengths of the merged reads

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `3_len_graphs/`
			- `<sample_name>_sequences_sizes_after_merging.png`:
			- `<sample_name>_sequences_sizes_after_merging_logscale.png`:	 

	</details>


- ### Step 4: Obtain sequences in fasta format
  
- ### Step 5: Obtain consensus sequences
  We have created an in-house Python (v3.6) script to calculate the consensus sequences by barcode
aiming to resolve PCR and RT errors.

First of all, we take the assembled reads and split it in order to obtain the sequence and the barcode
separately.
Secondly, we performed a cut-off process leaving out of the study barcodes which had less than the
lower cut-off value and more than the upper cut-off value sequences. With the rest of sequences, we
carried out a multiple alignment using [MAFFT](https://mafft.cbrc.jp/alignment/software/) software (Katoh et al., 2005). An absolute frequency matrix is
obtained with scores per base position in the consensus sequence. The resulted files from this
program execution are:

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `5_consensus/`
			- `<sample_name>_consensus.fna`: a FASTA file with the barcode, number of sequences and consensus sequence in each case
			- `<sample_name>_consensus.xls`: a summary of the results, showing the barcode, consensus and sequences per barcode	 
			- `<sample_name>_consensus.prf`: the matrix with specific score per position in the consensus sequence
			- `<sample_name>_discarded.txt`: a plain text file with sequences’ barcodes that didn’t fit the cutoff 
			- `<sample_name>_consensus.png`: a scatter plot of the barcode distribution (without the cutoff step)
			- `<sample_name>_cutoff_consensus.png`: a scatter plot of the barcode distribution (applying the cutoff step)
        
	</details>
 
- ### Step 6: Map consensus sequences
  We used the [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) software package which employs the Burrows-Wheeler Alignment tool and later we sorted the alignments with [SAMTools](https://github.com/samtools/samtools) package (Li et al., 2009), which provides diverse utilities for manipulating alignments in SAM format.

	<details>
	<summary>Output files</summary>
		
	- `Results/`
	   - `5_consensus/`
			- `<sample_name>_consensus.fna`: a FASTA file with the barcode, number of sequences and consensus sequence in each case
			- `<sample_name>_consensus.xls`: a summary of the results, showing the barcode, consensus and sequences per barcode	 
			- `<sample_name>_consensus.prf`: the matrix with specific score per position in the consensus sequence
			- `<sample_name>_discarded.txt`: a plain text file with sequences’ barcodes that didn’t fit the cutoff 
			- `<sample_name>_consensus.png`: a scatter plot of the barcode distribution (without the cutoff step)
			- `<sample_name>_cutoff_consensus.png`: a scatter plot of the barcode distribution (applying the cutoff step)
        
	</details>

- ### Step 7: Variant calling
  Aiming to calculate the frequency in which a variant is present in the consensus sequences compared to the reference genome we used the genetic variant detector software [FreeBayes](https://github.com/ekg/freebayes) (Garrison et al., 2012).
  
- ### Step 8: VCF analysis
  We have created an in-house Python (v3.6) script to calculate error rate from the VCF file obatined in the previous step as well as graphs.

	<details>
	<summary>Output files</summary>
		
	- `Results/`
	   - `8_vcf_analysis/`
			- `<sample_name>_variants_distribution.png`:
			- `<sample_name>_indels_distribution.png`:
			- `<sample_name>_heatmap_snp_types.png`:
			- `<sample_name>.xlsx`:
			- `<sample_name>.csv`:
   
	</details>


- ### Step 9: Offspring search (optional)
  Creates graphs similar to the ones described in Zhou et al., 2015.  

	<details markdown="1">
	<summary>Output files</summary>

	- `Results/`
	   - `9_offsprings/`
			- `<sample_name>_differences.png`:
   
	</details>



## Test data results

This is a summary of the expected results obtained running the pipeline with the test data (as indicated in [Quick Start](#quick-start)).

- ### Step 1: Quality Control
`<sample_name>_1_fastqc*.html`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/QC_1.png?raw=true)

- ### Step 3: Graph lengths of the merged reads
`<sample_name>_sequences_sizes_after_merging_logscale.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433_sequences_sizes_after_merging_logscale.png?raw=true)

- ### Step 5: Obtain consensus sequences
`<sample_name>_consensus.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433_consensus.png?raw=true)

- ### Step 8: VCF analysis
`<sample_name>_variants_distribution.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433_variants_distribution.png?raw=true)
  
`<sample_name>_indels_distribution.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433_indels_distribution.png?raw=true)
  
`<sample_name>_heatmap_snp_types.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433heatmap_snp_types.png?raw=true)
  
`<sample_name>.xlsx`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433_error_rate.png?raw=true)

- ### Step 9: Offspring search
`<sample_name>_differences.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_sample_output_images/O3MQ178433_differences.png?raw=true)

## Extra content
All the python scripts used in the pipeline are also available as independent Google Colab notebooks:

- [sizes-graphs](https://colab.research.google.com/github/friburgo-moc/FidelityFinder/blob/main/ColabNotebooks/sizes_graphs.ipynb?hl=es): creates plots with the lengths of the merged sequences obtained.
- [consensus-construction](https://colab.research.google.com/github/friburgo-moc/FidelityFinder/blob/main/ColabNotebooks/consensus_construction.ipynb?hl=es): finds barcodes and builds consensus sequences
- [vcf-analyzer](https://colab.research.google.com/github/friburgo-moc/FidelityFinder/blob/main/ColabNotebooks/vcf_analyzer.ipynb?hl=es): finds variants in VCF files. Creates an excel file with different data (table with variants, total number of variants, mutation rate...) and graphs showing the distribution of variants in the reference sequence, the distribution of indels and a heatmap with the types of SNPs (if any).
- [offsprings-finder](https://colab.research.google.com/github/friburgo-moc/FidelityFinder/blob/main/ColabNotebooks/offsprings_finder.ipynb?hl=es): finds possible offspring barcodes (barcodes generated due to errors in PCR reactions/sequencing reaction).

Additionally, there are two extra Google Colab notebooks to further process the results of the pipeline:
- [barcode-analyzer](https://colab.research.google.com/github/friburgo-moc/FidelityFinder/blob/main/ColabNotebooks/barcode_analyzer.ipynb?hl=es): creates a plot with a profile of the nucleotides present in barcodes
- [hotspots-finder](https://colab.research.google.com/github/friburgo-moc/FidelityFinder/blob/main/ColabNotebooks/hotspots_finder.ipynb?hl=es): it can process mutations tables of several samples to analyze and plot total variants, indels and SNPs in common

## Creative Commons
[![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/cc_logo.png?raw=true)](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)  
**Attribution-NonCommercial-ShareAlike 4.0 International**  
**(CC BY-NC-SA 4.0)** 

FidelityFinder includes several third-party packages provided under other open source licenses, please check them for additional details.

## Citation  
We politely request that this work be cited as:  

## References
Garrison E, G. Marth (2012) Haplotype-based variant detection from short-read sequencing.
arXiv:1207.3907 [q-bio.GN]

Katoh K., K. Kuma, H. Toh, T. Miyata (2005) MAFFT version 5: improvement in accuracy of
multiple sequence alignment. Nucleic Acids Research 33(2): 511–8

Li H., B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer, G. Marth, G. Abecasis, R.
Durbin, 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map
(SAM) format and SAMtools. Bioinformatics 25(16): 2078–2079

Zhang J., K. Kobert, T. Flouri, A. Stamatakis (2014) PEAR: A fast and accurate Illumina
Paired-End reAd merger. Bioinformatics 30(5): 614–620.

Zhou S., C. Jones, P. Mieczkowski, R. Swanstrom (2015) Primer ID Validates Template
Sampling Depth and Greatly Reduces the Error Rate of Next-Generation Sequencing of HIV-1
Genomic RNA Populations. Journal of Virology 89 (16) 8540-8555.
