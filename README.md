
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
    <a href="https://dl.circleci.com/status-badge/redirect/circleci/otEhDT8pGaAEyjBGSiC6d/5x57MPpNkWytxdZmFhA4Yb/tree/main">
        <img src="https://dl.circleci.com/status-badge/img/circleci/otEhDT8pGaAEyjBGSiC6d/5x57MPpNkWytxdZmFhA4Yb/tree/main.svg?style=shield&circle-token=7c1136d7caa0d43d54a9e8041af0a703d904348e" /></a>
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
   * [Step 4: Length and Quality Filtering and FASTQ to FASTA conversion](#step-4-length-and-quality-filtering-and-fastq-to-fasta-conversion)
   * [Step 5: Obtain Consensus Sequences](#step-5-obtain-consensus-sequences)
   * [Step 6: Map Consensus Sequences](#step-6-map-consensus-sequences)
   * [Step 7: Variant Calling](#step-7-variant-calling)
   * [Step 8: VCF Analysis](#step-8-vcf-analysis)
   * [Step 9: Offspring Search](#step-9-offspring-search)
6. [Test data results](#test-data-results)
7. [Extra content](#extra-content)
8. [Creative Commons](#creative-commons)
9. [Citation](#citation)
10. [Developers](#developers)
11. [References](#references)

## Introduction
**FidelityFinder** is a bioinformatics analysis pipeline to determine the fidelity of DNA synthesis by reverse transcriptases (RTs) from sequences obtained by Single Strand Consensus Sequencing (SSCS). 

Fidelity determination is based on the use of Unique Molecular Identifiers (UMIs) or barcodes. Each cDNA obtained by RTs is tagged with a barcode, so each cDNA molecule has a unique identity. Then, these cDNAs are amplified by PCR, and adapter sequences are added to generate a library (see [Libray preparation](#library-preparation) for a more detailed explanation). The libraries generated are then sequenced by NGS and **FidelityFinder** is able to evaluate the fidelity of the RT used: the pipeline discards PCR and NGS errors thanks to the construction of consensus sequences that share the same barcode sequence, and it obtains an error rate comparing the consensus sequences obtained and their reference sequence. This error rate is the combination of transcription and reverse transcription errors (see [Pipeline overview](#pipeline-overview) for a more detailed explanation of the pipeline). The higher the fidelity of the RT used, the lower the error rate.  

Single Strand Consensus Sequencing method to determine the fidelity of reverse transcriptases:

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Primers_IDs_method.PNG?raw=true)



## Quick Start 
1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=23.04.2)
2. Install [Docker](https://docs.docker.com/get-docker/) or [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html#regular-installation) 
3. Run the pipeline with the test data provided:
   ```console
   nextflow run https://github.com/RTlabCBM/FidelityFinder -profile <conda/docker>,test -r main
   ```
   > - Replace `<conda/docker>` for `conda` or `docker` depending on your election at step 2  
   > - You can find the expected results of the test data analysis in [Test data results](#test-data-results)
   
4. Run the pipeline with your own data:
   ```console
   nextflow run https://github.com/RTlabCBM/FidelityFinder -c <config_file> -profile <conda/docker>,<your_profile> -r main
   ```
   > - `<config_file>` must contain your input parameters for your own analysis. See [Input parameters](#input-parameters) section for more details.
   > - Replace `<conda/docker>` for `conda` or `docker` depending on your election at step 2
   > - Replace `<your_profile>` for the name of the profile with your own input parameters in the `config_file`
   > - If you are using conda, it is recommended to use the [`conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html#how-it-works) setting to store the environments in a central location for future pipeline runs.
   > - If you want to generate an [`Execution report`](https://www.nextflow.io/docs/latest/tracing.html#execution-report), add the following command line option: `-with-report [file name]`


## Input parameters

Input parameters must be provided in your own config file. You have to add your own profile to the `nextflow.config` file of this repository. It must contain the following structure:

  ```console
profiles {
	conda {
		conda.enabled = true
		docker.enabled = false
	}	
	docker {
		conda.enabled = false
		docker.enabled = true
		process.container = 'rtlabcbm/fidelityfinder_img2:latest' 		
	}
	your_profile  {
		params.seq_folder_path = "${baseDir}/raw_reads/"
		params.ref_seq_path = "${baseDir}/reference_sequences/reference.fasta"
		params.insert_length = "536"
		params.length_tolerance = "20"
		params.minimum_quality_score = "20"
		params.minimum_percent_quality_bases = "90"
		params.cutoff = "2"
		params.threshold = "50"
		params.min_pos = "23"
		params.max_pos = "492"
		params.cpus = "4"
 	}
}
  ```

- ### **params.seq_folder_path**
path of the folder where you have the raw sequences obtained by NGS. Files must have the following names: **<sample_name>_1.fastq** and **<sample_name>_2.fastq** for the forward and reverse sequences, respectively. You can also have them compressed: **<sample_name>_1.fastq.gz** and **<sample_name>_2.fastq.gz**. Example:
```console
params.seq_folder_path = "${baseDir}/raw_reads/"
```

- ### **params.ref_seq_path**
path of the file with the reference sequence, i.e., the sequence of the insert of the library (without mutations). The barcode/s sequence/s must be indicated with as many "N" as nucleotides the barcode has. It should be written in the same 5'-->3' sense of the forward reads (or Reads 1s) of the NGS result. Example:
```console
params.ref_seq_path = "${baseDir}/reference_sequences/reference.fasta"
```

And this would be an example of the content of the fasta file with an insert sequence that has a barcode of 14 nucleotides:

```console
>insert_reference_sequence
CTTCCTACAAGGGAATTGGAGGTGGAATGGATGGCCCAAAAGTTAAACNNNNNNNNNNNNNNACCTT
```

- ### **params.insert_length**
length (integer number of nucleotides) of the library insert that has been sequenced, i.e., total length of the library except for the adapter sequences. Example:
```console
params.insert_length = "536"
```

- ### **params.length_tolerance**
length (integer number of nucleotides) for filtering merged reads. This parameter is important to filter merged reads according to their length: reads that differ by more than **params.length_tolerance** nucleotides from the indicated **params.insert_length** are not selected. For example, if params.insert_length is set to 536, and params.length_tolerance is set to 20, this means that merged reads with lengths between 516 and 556 nucleotides (536 ± 20) will be considered for further analysis. Example:
```console
params.insert_length = "20"
```

- ### **params.minimum_quality_score**
minimum quality score (integer number) required to keep a sequence (also known as Phred quality score). Sequences with quality scores below this threshold will be filtered out. Use a value of 0 to retain all sequences and avoid quality filtering. Example:
```console
params.cutoff = "10"
```

- ### **params.minimum_percent_quality_bases**
minimum percent of bases (integer number) that must have the specified quality score (**params.minimum_quality_score**). If the percentage of bases in a sequence with a quality score equal to or greater than **params.minimum_quality_score** is below this threshold, the sequence will be filtered out. Use a value of 0 to retain all sequences and avoid quality filtering. Example:
```console
params.cutoff = "90"
```

- ### **params.cutoff**
this cutoff (integer number) is used to discard reads during the consensus construction: if the number of reads that share a barcode is equal to or lower than the cutoff, these reads will not be used. A minimum of 3 reads with the same barcode is needed to build an appropriate consensus sequence, so the minimum cutoff value should be 2. Example:
```console
params.cutoff = "2"
```

- ### **params.threshold**
threshold percentage (integer number in the range 0-100) is used to determine the consensus sequences of several reads with the same barcode. Sequences that share a barcode are aligned and, for each position of the alignment, the proportion of each nucleotide (or deletion) is calculated. If the proportion of one nucleotide (or deletion) is equal to or higher than the threshold, the nucleotide (or deletion) is added to the consensus sequence of the aligned reads, otherwise an "N" will be incorporated and not taken into account as a position with an error.

The threshold parameter admits values between 0 and 100. For example, if the threshold is 90 and there are 10 reads with the same barcode, for each position, the consensus sequence is built using the nucleotides (or deletions) present in at least 9 of the reads, otherwise an "N" will be incorporated. Example:

```console
params.threshold = "90"
```

A threshold of 100 would be more strict, each position must have the same nucleotide (or deletion) in all the aligned reads; while a threshold of 0 would mean that the consensus sequence is built using the nucleotides (or deletions) that are present in the majority of the aligned reads, for each position. If for a given position there is no majority nucleotide (or deletion), e.g. in 50% of the reads there is a "T" and in the other 50% a "C",  an "N" is always added to the consensus sequence, regardless of the chosen threshold.

- ### **params.min_pos**
the first position (integer number) of the reference sequence used to quantify mutations during the VCF analysis step. This parameter is useful to not consider the beginning of the library insert in case it contains a sequence that does not come directly from the cDNA synthesized during the reverse transcription. For example, if the first 15 nucleotides of your insert are a primer binding sequence during the library preparation and/or contain a barcode, params.min_pos value should be 16. Example:
```console
params.min_pos = "16"
```

- ### **params.max_pos**
last position (integer number) of the reference sequence used to quantify mutations during the VCF analysis step. This parameter is useful to not consider the end of the library insert in case it contains a sequence that does not come directly from the cDNA synthesized during the reverse transcription. For example, if the last 15 nucleotides (of an insert with a total length of 100 nucleotides) are a primer binding sequence during the library preparation and/or contain a barcode, params.max_pos value should be 84. Example:
```console
params.max_pos = "84"
```

- ### **params.cpus**
number of CPU cores (integer number) available to run the nextflow pipeline. This parameter specifies the computational resources the program can utilize for parallel processing. Set to 1 if you only want to use a single CPU core for its computations.
```console
params.cpus = "1"
```

Once you have created your config file with your own parameters, you have to run the pipeline and specify the path of the config_file with `-c <config_file>`, and the name of the profile with your specific parameters: `-profile my_profile`.

## Library preparation

One way to prepare the libraries is as follows: 
1. **Reverse transcription** step using a primer with an overhang that contains a barcode sequence
2. **PCR** step to amplify obtained cDNAs and incorporate the adapter sequences
3. **Next Generation Sequencing**: short-read sequencing technology in paired-end mode
4. **Analysis** of the fastq files (forward and reverse reads, or Reads 1 and Reads 2) with [FidelityFinder](#pipeline-overview )

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Library_preparation_workflow.PNG?raw=true)

It is possible to vary the way the libraries are prepared. For example, by adding several barcodes, adding barcodes by ligation, using intermediate purification steps... 

## Pipeline overview 

![Workflow](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/Workflow.PNG?raw=true)

The pipeline is built using Nextflow. Processing steps:

- ### Step 1: Quality Control
  Quality analyses are performed over reads (forward and reverse) using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) software. It exports the results to an HTML-based permanent report.

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `1_quality_control/`
    		- `<sample_name>_1_fastqc*.html`: html file
      		- `<sample_name>_1_fastqc*.zip`: zip file	 
    		- `<sample_name>_2_fastqc*.html`: html file
      		- `<sample_name>_2_fastqc*.zip`: zip file

	</details>

  
- ### Step 2: Joining of paired reads
  Paired reads are joined into one single sequence. This is performed using [PEAR](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (Paired-End reAd mergeR) (Zhang *et al.*, 2014), a fast and accurate paired-end read merger. PEAR evaluates all possible paired-end read overlaps without requiring the target fragment size as input. In addition, it implements a statistical test for minimizing false-positive results. PEAR outputs four files: a file containing the assembled reads, two files containing the forward and reverse unassembled reads, and a file containing the discarded reads.
	
	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `2_merged_reads/`
			- `<sample_name>.assembled.fastq`: file containing the assembled reads by PEAR
			- `<sample_name>.unassembled.forward.fastq`: file containing forward unassembled reads by PEAR
			- `<sample_name>.unassembled.reverse.fastq`: file containing reverse unassembled reads by PEAR
			- `<sample_name>.discarded.fast`: file containing discarded reads by PEAR
			- `<sample_name>_pear.log`: log file of the PEAR program
	</details>

- ### Step 3: Graph lengths of the merged reads
  An in-house Python (v3.6) script is used to graph the lengths of the assembled reads of the [Step 2](#step-2-joining-of-paired-reads-and-filtering-by-length). It uses the `<sample_name>._lengths.txt` file as input.

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `3_reads_lengths_graphs/`
    			- `<sample_name>._lengths.txt`: file with info about the assembled reads lengths. The first column is the length (in nucleotides) and the second column indicates the number of reads with the specified length
			- `<sample_name>_lengths_distribution.png`: graph with the lengths of the merged reads of the [Step 2](#step-2-joining-of-paired-reads-and-filtering-by-length)
			- `<sample_name>_lengths_distribution_logscale.png`: graph with the lengths of the merged reads of the [Step 2](#step-2-joining-of-paired-reads-and-filtering-by-length) using a log scale

	</details>

- ### Step 4: Length and quality filtering and FASTQ to FASTA conversion
  Reads from the `<sample_name>.assembled_filtered.fastq` file of the [Step 2](#step-2-joining-of-paired-reads) are filtered.
  Length filtering:
  	Reads that differ by more than **params.length_tolerance** nucleotides from the **params.insert_length** are not selected.
  Quality filtering:
  	Reads are filtered according to the quality values specified in **params.minimum_quality_score** and **params.minimum_percent_quality_bases**. [FASTQ Quality Filter](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage) tool is used. A quality control (see [Step 1](#step-1-quality-control)) is performed before and after quality filters. 
  FASTQ to FASTA conversion:
  After filtering, selected reads are then converted into FASTA format using the [FASTQ-to-FASTA](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_to_fasta_usage) tool. 

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `4_fasta_sequences_quality_filtered/`
			- `<sample_name>_fastq_quality_filter.log`: log file of the FASTQ Quality Filter program
			- `<sample_name>_fastq_to_fasta.log`: log file of the FASTQ-to-FASTA program  
			- `<sample_name>_lengths_and_qual_filtered.fasta`: fasta file with the selected reads after filtering by length and sequencing quality
        		- `<sample_name>_lengths_and_qual_filtered.fastq`: fastq file with the selected reads after filtering by length and sequencing quality
  			- `<sample_name>_lengths_filtered.fastq`: fastq file with the selected reads after filtering by length     
			- `<sample_name>_lengths_and_qual_filtered_fastqc.html`: html file with quality information after length and sequencing quality filtering
			- `<sample_name>_lengths_and_qual_filtered_fastqc.zip`: zip file with quality information after length and sequencing quality filtering
 			- `<sample_name>_lengths_filtered_fastqc.html`: html file with quality information after length filtering
			- `<sample_name>_lengths_filtered_fastqc.zip`: zip file with quality information after length filtering    
	</details>

  
- ### Step 5: Obtain consensus sequences
  An in-house Python (v3.6) script is used to calculate the consensus sequences by barcode aiming to resolve PCR and NGS errors. First, the sequences from the `<sample_name>_qual_filtered.fasta` file of the [Step 4](#step-4-quality-filtering-and-fastq-to-fasta-conversion) are aligned with respect to the reference sequence. The barcode sequence is identified by matching the nucleotides marked as "N" in the reference sequence. If the barcode identified has the same length as the barcode of the reference sequence, the read is selected. If there is more than one barcode in the reference sequence, the identified barcode will be the concatenation of them. 

  Secondly, the sequences that share the same barcode are grouped together. If the number of sequences with the same barcode is equal to or lower than the input cutoff value (**params.cutoff**), they are discarded. The selected reads sharing a barcode are then aligned using [MAFFT](https://mafft.cbrc.jp/alignment/software/) software (Katoh *et al.*, 2005), and a consensus sequence is constructed using the threshold indicated as input (**params.threshold**).

	Among the output files, a fasta file is generated with the obtained consensus sequences (used in the [next step](#step-6-map-consensus-sequences)) and a JSON file with the sequences of the identified barcodes and their frequencies (used in the [Step 9](#step-9-offspring-search) to search for offspring barcodes).

	<details markdown="1">
	<summary>Output files</summary>
		
	- `Results/`
	   - `5_consensus_sequences/`
     			- `<sample_name>_barcodes.json`: a JSON file with the sequences of the identified barcodes and their frequencies      
			- `<sample_name>_consensus.prf`: a matrix with specific score per position in the consensus sequence
			- `<sample_name>_consensus_construction.log`: log file of the consensus_construction.py program
			- `<sample_name>_discarded.txt`: a plain text file with sequences barcodes that did not fit the cutoff 
			- `<sample_name>_frequencies_distribution_graph.png`: a scatter plot of the barcode distribution
			- `<sample_name>_t<params.threshold>_consensus.fna`: a FASTA file with the consensus sequences identified, the header of each sequence contains the barcode and the number of sequences used for consensus construction
	</details>
 

 
- ### Step 6: Map consensus sequences
  Reads from the `<sample_name>_consensus.fna` file of the [Step 5](#step-5-obtain-consensus-sequences) are aligned to the reference sequence using the [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) software package, which employs the Burrows-Wheeler Alignment tool. Afterwards, alignments are sorted with [SAMTools](https://github.com/samtools/samtools) package (Li *et al.*, 2009), which provides diverse utilities for manipulating alignments in SAM format.

	<details>
	<summary>Output files</summary>
		
	- `Results/`
	   - `6_mapped_consensus_sequences/`
			- `<sample_name>_t<params.threshold>_align.sam`: Sequence Alignment/Map (SAM) file of the consensus sequences
    			- `<sample_name>_t<params.threshold>_align_sort`: Binary Alignment Map file of the consensus sequences

        
	</details>

- ### Step 7: Variant calling
  The variant detector software [FreeBayes](https://github.com/ekg/freebayes) (Garrison *et al.*, 2012) is used to find the variants present in the consensus sequences. It uses the `<sample_name>.bam` file of the [Step 6](#step-6-map-consensus-sequences) as input.

  	<details>
	<summary>Output files</summary>
		
	- `Results/`
	   - `7_vcf_file/`
			- `<sample_name>_t<params.threshold>.vcf`: Variant Call Format (VCF) file

	</details>
  
- ### Step 8: VCF analysis
  An in-house Python (v3.6) script is used to analyze the variants information of the `<sample_name>.vcf` file of the [Step 7](#step-7-variant-calling). The script creates a report (an Excel file) with different data (table with variants, total number of variants, error rate...) and creates graphs showing the distribution of variants in the reference sequence, the distribution of indels, and a heatmap with the types of SNPs (if any).
  
	<details>
	<summary>Output files</summary>
		
	- `Results/`
	   - `8_vcf_analysis/`
			- `<sample_name>_t<params.threshold>_variants_distribution.png`: a graph showing the number of variants detected in each position of the reference sequence
			- `<sample_name>_t<params.threshold>_indels_distribution.png`: a graph showing the number of indels detected in each position of the reference sequence
			- `<sample_name>_t<params.threshold>_heatmap_snp_types.png`: a heatmap chart with the proportion of the different SNP types detected
			- `<sample_name>_t<params.threshold>.xlsx`: an Excel file with different sheets containing tables and statistics of the variants detected
			- `<sample_name>_t<params.threshold>.csv`: a csv file with a list the variants detected in each position of the reference sequence
   
	</details>


- ### Step 9: Offspring search
  An in-house Python (v3.6) script is used to identify possible offspring barcodes. It uses the `<sample_name>_barcodes.json` of the [Step 5](#step-5-obtain-consensus-sequences) as input. It follows a similar strategy to the one described in Zhou *et al.*, 2015. Two types of offspring barcodes are identified: **barcodes with 1 difference** with respect to other barcodes of higher frequency and **barcodes with 2 differences** with respect to other barcodes of higher frequency. This step can be slow if the number of barcodes is high. It can be omitted if a faster analysis is desired.

	<details markdown="1">
	<summary>Output files</summary>

	- `Results/`
	   - `9_offsprings/`
			- `<sample_name>_differences.png`: a graph showing the distribution of barcodes according to their frequency and the distribution of barcodes with 1 or 2 differences with respect to other barcodes of higher frequency.
			- `<sample_name>_percentage_differences.png`: a graph showing the distribution of barcodes using a percentages to show offspring barcodes
			- `<sample_name>_offsprings_finder.log`: log file of the offsprings_finder.py program
	</details>
 
>[!TIP]
>The results of the last step (Offspring search) may be useful to assess whether the cutoff chosen in the analysis was appropriate, or to consider repeating the analysis with a new cutoff value to avoid taking into account offspring barcodes.

## Test data results

This is a summary of the expected results obtained running the pipeline with the test data (as indicated in [Quick Start](#quick-start)).

- ### Step 3: Graph lengths of the merged reads
`test_frequencies_distribution_graph.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_output_images/test_lengths_distribution_logscale.png?raw=true)

- ### Step 5: Obtain consensus sequences
`<sample_name>_consensus.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_output_images/test_frequencies_distribution_graph.png?raw=true)

- ### Step 8: VCF analysis
`test_t50_variants_distribution.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_output_images/test_t50_variants_distribution.png?raw=true)
  
`test_t50heatmap_snp_types.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_output_images/test_t50heatmap_snp_types.png?raw=true)
  
`test_t50.xlsx`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_output_images/test_t50.xlsx.png?raw=true)

- ### Step 9: Offspring search
`test_percentage_differences.png`  
![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/test_output_images/test_percentage_differences.png?raw=true)

## Extra content
All the Python scripts used in the pipeline are also accessible as independent Jupyter Notebooks (available at [FideltyFinderJupyter](https://github.com/RTlabCBM/FidelityFinderJupyter)). They can be opened with the following Google Colab links:

- [lengths_distribution_plotting](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/lengths_distribution_plotting.ipynb): creates plots with the lengths of the merged sequences obtained.
- [consensus_construction](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/consensus_construction.ipynb): finds barcodes and builds consensus sequences
- [vcf_analyzer](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/vcf_analyzer.ipynb): finds variants in VCF files. Creates an excel file with different data (table with variants, total number of variants, mutation rate...) and graphs showing the distribution of variants in the reference sequence, the distribution of indels and a heatmap with the types of SNPs (if any).
- [offsprings_finder](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/offsprings_finder.ipynb): finds possible offspring barcodes (barcodes generated due to errors in PCR reactions/sequencing reaction).

Additionally, there are two extra Google Colab notebooks to further process the results of the pipeline and simulate the pipeline:
- [barcode_analyzer](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/barcode_analyzer.ipynb): creates a plot with a profile of the nucleotides present in barcodes
- [hotspots_finder](https://colab.research.google.com/github/RTlabCBM/FidelityFinderJupyter/blob/main/JupyterNotebooks/hotspots_finder.ipynb): it can process mutations tables of several samples to analyze and plot total variants, indels and SNPs in common

## Creative Commons
[![image](https://github.com/RTlabCBM/FidelityFinder/blob/main/docs/images/cc_logo.png?raw=true)](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)  
**Attribution-NonCommercial-ShareAlike 4.0 International**  
**(CC BY-NC-SA 4.0)** 

FidelityFinder includes several third-party packages provided under other open source licenses, please check them for additional details.

## Citation  
We politely request that this work be cited as:

Martínez del Río, J., Frutos-Beltrán, E., Sebastián-Martín, A., Lasala, F., Yasukawa, K., Delgado, R., & Arias, L. M. (2024). HIV-1 reverse transcriptase error rates and transcriptional thresholds based on single-strand consensus sequencing of target RNA derived from in vitro-transcription and HIV-infected cells. Journal of Molecular Biology, 168815. [https://doi.org/10.1016/j.jmb.2024.168815](https://doi.org/10.1016/j.jmb.2024.168815)

## Developers
### Main developer
- Javier Martínez del Río
  
### Collaborators
- Estrella Frutos-Beltrán

## References
Garrison E, G. Marth (2012) Haplotype-based variant detection from short-read sequencing.
arXiv:1207.3907 [q-bio.GN].

Katoh K., K. Kuma, H. Toh, T. Miyata (2005) MAFFT version 5: improvement in accuracy of
multiple sequence alignment. Nucleic Acids Research 33(2): 511–8.

Li H., B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer, G. Marth, G. Abecasis, R.
Durbin, 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map
(SAM) format and SAMtools. Bioinformatics 25(16): 2078–2079.

Zhang J., K. Kobert, T. Flouri, A. Stamatakis (2014) PEAR: A fast and accurate Illumina
Paired-End reAd merger. Bioinformatics 30(5): 614–620.

Zhou S., C. Jones, P. Mieczkowski, R. Swanstrom (2015) Primer ID Validates Template
Sampling Depth and Greatly Reduces the Error Rate of Next-Generation Sequencing of HIV-1
Genomic RNA Populations. Journal of Virology 89 (16) 8540-8555.
