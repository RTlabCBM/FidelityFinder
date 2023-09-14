
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
**FidelityFinder** is a bioinformatics analysis pipeline to determine the fidelity of transcriptases and the fidelity of DNA synthesis by reverse transcriptases (RTs) from sequences obtained by Next Generation Sequencing (NGS). 

Fidelity determination is based on the "Primers IDs method". Each cDNA obtained by RTs is tagged with a barcode, so each cDNA molecule has a unique identity. Then, these cDNAs are amplified by PCR and adaptor sequences are added to generate a library (see [Libray preparation](#library-preparation) for a more detailed explanation). The libraries generated are sequenced then by NGS and **FidelityFinder** is able to evaluate the fidelity of the RT used: the pipeline is able to discard PCR and NGS errors thanks to the construction of consensus sequences that share the same barcode sequence, so it obtains and error rate. This error rate is the combination of transcription and reverse transcription errors (see [Pipeline overview](#pipeline-overview) for a more detailed explanation of the pipeline). The higher the fidelity of the RT used, the lower the error rate.  

Primers IDs method to determine the fidelity of reverse transcriptases:

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Primers_IDs_method.PNG?raw=true)

## Quick Start 
1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=23.04.2)
2. Install [Docker](https://docs.docker.com/get-docker/) or [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html#regular-installation) 
3. Run the pipeline with the test data provided:
   ```console
   nextflow run https://github.com/RTlabCBM/FidelityFinder -profile <conda/docker>,test_profile -r main
   ```
   > - Replace `<conda/docker>` for `conda` or `docker` depending on your election at step 2  
   > - You can find the expected results of the test data analysis in [Test data results](#test-data-results)
   
4. Run the pipeline with your own data:
   ```console
   nextflow run https://github.com/RTlabCBM/FidelityFinder -c <config_file> -<conda/docker>,<your_profile> -r main
   ```
   > - `<config_file>` must contain your input parameters for your own analysis. See [Input parameters](#input-parameters) section for more details.
   > - Replace `<conda/docker>` for `conda` or `docker` depending on your election at step 2
   > - Replace `<your_profile>` for the name of the profile with your own input parameters in the `config_file`
   > - If you are using conda, it is recommended to use the [`conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html#how-it-works) setting to store the environments in a central location for future pipeline runs.
   > - If you want to generate an [`Execution report`](https://www.nextflow.io/docs/latest/tracing.html#execution-report), add the following command line option: `-with-report [file name]`


## Input parameters

Input parameters must be provided in your own config file. You have to add your own profile to the `nextflow.config` file of this repository. It must contain the following structure:

  ```console
 your_profile  {
	params.seq_folder_path = "/folder/with/raw/sequences"
	params.ref_seq_path = "/file/with/reference/sequence"
	params.insert_length = "length of the library insert"
	params.cutoff = "cutoff_value_to_discard_low_frequency_barcodes"
	params.bc_size = "nucleotides_of_the_barcode"
	params.threshold = "threshold_used_to_construct_consensus_sequences"
	params.min_pos = "first position of the reference sequence used to quantify mutations"
	params.max_pos = "last position of the reference sequence used to quantify mutations"
 }
  ```

- ### **params.seq_folder_path**
path of the folder where you have the raw sequences obtained by NGS. Files must have the following names: **<sample_name>_1.fastq** and **<sample_name>_2.fastq** for the forward and reverse sequences, respectively. You can also have them compressed: **<sample_name>_1.fastq.gz** and **<sample_name>_2.fastq.gz**. Example:
```console
params.seq_folder_path = "${baseDir}/RawReads/"
```

- ### **params.ref_seq_path**
path of the file with the reference sequence, i.e., the sequence of the insert of the library (without mutations). The barcode/s sequence/s must be indicated with as many "N" as nucleotides the barcode has. It should be written in the same 5'-->3' sense of the forward read (or read 1) of the NGS result. Example:
```console
params.ref_seq_path = "${baseDir}/Sequences/insert_reference_sequence.fasta"
```

<p> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; And this would be an example of the content of the fasta file with an insert sequence that has a barcode of 14 nucleotides: </p>   

```console
>insert_reference_sequence
CTTCCTACAAGGGAATTGGAGGTGGAATGGATGGCCCAAAAGTTAAACNNNNNNNNNNNNNNACCTT
```

- ### **params.insert_length (HARÉ QUE SE DETERMINE AUTOMÁTICAMENTE CONTANDO EN NÚMERO DE NUCLEÓTIDOS PROPORCIONADO EN ref_seq_path)**
length of the library insert that has been sequenced, i.e. total length of the library except for the adaptors. This parameter is important to filter merged reads according to their length: reads that differ by more than 20 nucleotides from the indicated length are not selected. 

- ### **params.cutoff**
this cutoff is used to discard reads during the consensus construction: if the number of reads that share a barcode is equal to or lower than the cutoff, these reads will not be used. A minimum of 3 reads with the same barcode is needed to build an appropriate consensus sequence, so the minimum cutoff value should be 2. Example:
```console
params.cutoff = "2"
```

- ### **params.bc_size (HARÉ QUE SE DETERMINE AUTOMÁTICAMENTE CONTANDO EN NÚMERO DE Ns PROPORCIONADO EN ref_seq_path)**

- ### **params.threshold**
this threshold value is used to determine the consensus sequences of several reads with the same barcode. Sequences that share a barcode are aligned, and for each position of the alignment it is calculated the proportion of each nucleotide (or deletion). If the proportion of one nucleotide (or deletion) is equal to or higher than the threshold, the nucleotide (or deletion) is added to the consensus sequence of the aligned reads, otherwise an "N" will be incorporated and not taken into account as a position with an error.

The threshold parameter admits values between 0 and 1. For example, if the threshold is 0.9 and there are 10 reads with the same barcode, the consensus sequence is built using the nucleotides (or deletions) present in at least 9 of the reads for each position, otherwise an "N" will be incorporated. A threshold of 1 would be more strict, each position must have the same nucleotide (or deletion) in all the aligned reads; while a threshold of 0 would mean that the consensus sequence is built using the nucleotides (or deletion) that are present in the majority of the aligned reads, for each position. If for a given position there is no majority nucleotide (or deletion), e.g. in 50% of the reads there is a "T" and in the other 50% a "C",  an "N" is always added to the consensus sequence, regardless of the chosen threshold. Example:

```console
params.threshold = "0.9"
```

- ### **params.min_pos**
first position of the reference sequence used to quantify mutations during the VCF analysis step. This parameter is useful to not consider the beginning of the library insert in case it contains a sequence that does not come directly from the cDNA synthesized during the reverse transcription. For example, if the first 15 nucleotides of your insert are a primer binding sequence during the library preparation and/or contains a barcode, params.min_pos value should be 16. Example:
```console
params.min_pos = "16"
```

- ### **params.max_pos**
last position of the reference sequence used to quantify mutations during the VCF analysis step. This parameter is useful to not consider the end of the library insert in case it contains a sequence that does not come directly from the cDNA synthesized during the reverse transcription. For example, if the last 15 nucleotides (of an insert with a total length of 100 nucleotides) are a primer binding sequence during the library preparation and/or contains a barcode, params.max_pos value should be 84. Example:
```console
params.max_pos = "84"
```

Once you have created your config file with your own parameters, you have to run the pipeline and specify the path of the config_file with -c <config_file> and the name of the profile with your specific parameters -profile my_profile

## Library preparation

One way to prepare the libraries is as follows: 
1. **Reverse transcription** step using a primer with an overhang that contains a barcode sequence
2. **PCR** step to incorporate the adapter sequences
3. **Next Generation Sequencing**: short-read sequencing technology in paired-end mode
4. **Analysis** of the fastq files (forward and reverse reads, or read 1 and read 2) with [FidelityFinder](#pipeline-overview )

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Library_preparation_workflow.PNG?raw=true)

It is possible to vary the way the libraries are prepared. For example, by adding several barcodes, adding barcodes by ligation, using intermediate purification steps... 

## Pipeline overview 

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Workflow.PNG?raw=true)

The pipeline is built using Nextflow. Processing steps:

- ### Step 1: Quality Control
  Quality analyses are performed over reads using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) software. It exports the results to an HTML-based permanent report.

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `1_fastqc_results/`
    		- `<sample_name>_1_fastqc*.html`: html file
      		- `<sample_name>_1_fastqc*.zip`: zip file	 
    		- `<sample_name>_2_fastqc*.html`: html file
      		- `<sample_name>_2_fastqc*.zip`: zip file

	</details>

  
- ### Step 2: Joining of paired reads and filtering by length
  Paired reads are joined into one single sequence. This is performed using [PEAR](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (Paired-End reAd mergeR) (J. Zhang et al., 2014), a fast and accurate paired-end read merger. PEAR evaluates all possible paired-end read overlaps without requiring the target fragment size as input. In addition, it implements a statistical test for minimizing false-positive results. PEAR outputs four files: A file containing the assembled reads (`<sample_name>.assembled.fastq`), two files containing the forward and reverse unassembled reads (`<sample_name>.unassembled.forward.fastq` and `<sample_name>.unassembled.reverse.fastq`), and a file containing the discarded reads (`<sample_name>.discarded.fast`). Merged reads are then filtered according to their length (`<sample_name>.assembled_filtered.fastq`): reads that differ by more than 20 nucleotides from the reference insert length are not selected. 
	
	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `3_len_graphs/`
			- `<sample_name>.assembled.fastq`:
			- `<sample_name>.unassembled.forward.fastq`:
			- `<sample_name>.unassembled.reverse.fastq`:
			-  <sample_name>.discarded.fast`:
 			- `<sample_name>.assembled_filtered.fastq`:
			- `<sample_name>._lengths.txt`:
	</details>

				
	 If you want to modify the ±20 nucleotides filter, you have to change the "20" numbers found in the following line of `main.nf`:

```console
awk 'BEGIN {FS = "\\t" ; OFS = "\\n"} {header = \$0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= ("${params.insert_length}"-20) && length(seq) <= ("${params.insert_length}"+20)) {print header, seq, qheader, qseq}}' "${sampleId}.assembled.fastq" > "${sampleId}.assembled_filtered.fastq"
```

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
  An in-house Python (v3.6) script is used to calculate the consensus sequences by barcode aiming to resolve PCR and RT errors.

  First of all, we take the assembled reads and split it in order to obtain the sequence and the barcode separately.

  Secondly, we performed a cut-off process leaving out of the study barcodes which had less than the lower cut-off value and more than the upper cut-off value sequences. With the rest of the sequences, we carried out a multiple alignment using [MAFFT](https://mafft.cbrc.jp/alignment/software/) software (Katoh et al., 2005). An absolute frequency matrix is obtained with scores per base position in the consensus sequence. The resulting files from this program execution are:
  
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
	   - `6_map_consensus/`
			- `<sample_name>.sam`: a sam file
    			- `<sample_name>.bam`: a bam file

        
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
All the Python scripts used in the pipeline are also accessible as independent Jupyter Notebooks (available at [FideltyFinderJupyter](https://github.com/RTlabCBM/FidelityFinderJupyter)). They can be opened with the following Google Colab links:

- [sizes-graphs](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/sizes_graphs.ipynb?hl=es): creates plots with the lengths of the merged sequences obtained.
- [consensus-construction](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/consensus_construction.ipynb?hl=es): finds barcodes and builds consensus sequences
- [vcf-analyzer](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/vcf_analyzer.ipynb?hl=es): finds variants in VCF files. Creates an excel file with different data (table with variants, total number of variants, mutation rate...) and graphs showing the distribution of variants in the reference sequence, the distribution of indels and a heatmap with the types of SNPs (if any).
- [offsprings-finder](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/offsprings_finder.ipynb?hl=es): finds possible offspring barcodes (barcodes generated due to errors in PCR reactions/sequencing reaction).

Additionally, there are two extra Google Colab notebooks to further process the results of the pipeline:
- [barcode-analyzer](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/barcode_analyzer.ipynb?hl=es): creates a plot with a profile of the nucleotides present in barcodes
- [hotspots-finder](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/hotspots_finder.ipynb?hl=es): it can process mutations tables of several samples to analyze and plot total variants, indels and SNPs in common

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
