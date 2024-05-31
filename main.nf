#!/usr/bin/env nextflow


/* 

	############################################################
	###################### FidelityFinder ######################
	############################################################
	
	DESCRIPTION
	FidelityFinder is a bioinformatics analysis pipeline to determine the fidelity of DNA synthesis by 
	reverse transcriptases (RTs) from sequences obtained by Next Generation Sequencing (NGS).
	
	For additional details, please refer to the GitHub repository: https://github.com/RTlabCBM/FidelityFinder
	
	AUTHOR
	Javier Martínez del Río (javier.martinez@cbm.csic.es; javier.mardelrio@gmail.com)
	
	LICENSE
	Attribution-NonCommercial-ShareAlike 4.0 International
	(CC BY-NC-SA 4.0) (https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)

	FidelityFinder includes several third-party packages provided under other open source licenses, please check them for additional details.
	
*/


nextflow.enable.dsl=2


// Step 1: Quality Control
process quality_control {

	publishDir 'Results/1_quality_control', mode: 'copy'

	conda "${baseDir}/conda/processing_env.yml"
	
    input:
    path reads_file
	
	output:
	path "*"
	
    script:
    """
    fastqc ${reads_file}
    """
}


// Step 2: Joining of paired reads
process reads_merging {

    publishDir 'Results/2_merged_reads', mode: 'copy'
	
	conda "${baseDir}/conda/processing_env.yml"
    
	input:
	tuple val(sampleID), path(reads_files)
	
	output:   
    path "${sampleID}.unassembled.forward.fastq"
    path "${sampleID}.unassembled.reverse.fastq"
    path "${sampleID}.discarded.fastq"
	path "${sampleID}.assembled.fastq"
    path "${sampleID}_pear.log"
	tuple val(sampleID), path("${sampleID}.assembled.fastq"), emit: assembled_reads_ch

	
	script:
    """
	# Assemble paired-end reads using PEAR
		pear -f ${reads_files[0]} -r ${reads_files[1]} -j ${task.cpus} -o "${sampleID}" | tee "${sampleID}_pear.log"  
	"""
}


// Step 3: Graph lengths of the merged reads
process graph_reads_lengths {

	publishDir 'Results/3_reads_lengths_graphs', mode: 'copy'
	
	conda "${baseDir}/conda/python_env.yml"
	
	input:
	tuple val(sampleID), path(assembled_reads_file)
	
    output:
	path "${sampleID}_lengths.txt"
	path "${sampleID}_lengths_distribution.png"
	path "${sampleID}_lengths_distribution_logscale.png"
	
	script:
    """
	# Analyze assembled sequence lengths and save to txt file
		awk 'NR % 4 == 2 {dict_lengths[length(\$0)]++} END {for (l in dict_lengths) {print l, dict_lengths[l]}}' "${assembled_reads_file}" | sort -n > "${sampleID}_lengths.txt"
	# Create bar plots depicting the distribution of sequence lengths	
		python ${baseDir}/bin/lengths_distribution_plotting.py --input_file "${sampleID}_lengths.txt" --output_prefix "${sampleID}"
	"""
}


// Step 4: Length and quality Filtering | FASTQ to FASTA conversion
process quality_filtering_and_fastq_to_fasta {

	publishDir 'Results/4_fasta_sequences_quality_filtered', mode: 'copy'
	
	conda "${baseDir}/conda/processing_env.yml"
	
	input:
	tuple val(sampleID), path(assembled_reads_fastq)
    val params.insert_length
	val params.length_tolerance
	val params.minimum_quality_score
	val params.minimum_percent_quality_bases
	
    output:
	path "${sampleID}_lengths_filtered.fastq"
	path "${sampleID}_lengths_and_qual_filtered.fastq"
	path "${sampleID}_fastq_quality_filter.log"
	tuple val(sampleID), path("${sampleID}_lengths_and_qual_filtered.fasta"), emit: fasta_file_ch
	path "${sampleID}_fastq_to_fasta.log"
	path "${sampleID}_lengths_filtered_fastqc.html"
	path "${sampleID}_lengths_filtered_fastqc.zip"
	path "${sampleID}_lengths_and_qual_filtered_fastqc.html"
	path "${sampleID}_lengths_and_qual_filtered_fastqc.zip"
	
	script:
    """	
	# Filter assembled sequences by length (params.insert_length +- params.length_tolerance nt) and by quality
		awk 'BEGIN {FS = "\\t" ; OFS = "\\n"} {header = \$0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= ("${params.insert_length}"-"${params.length_tolerance}") && length(seq) <= ("${params.insert_length}"+"${params.length_tolerance}")) {print header, seq, qheader, qseq}}' "${assembled_reads_fastq}" > "${sampleID}_lengths_filtered.fastq"
		fastqc "${sampleID}_lengths_filtered.fastq"
		fastq_quality_filter -v -q "${params.minimum_quality_score}" -p "${params.minimum_percent_quality_bases}" -i "${sampleID}_lengths_filtered.fastq" -o "${sampleID}_lengths_and_qual_filtered.fastq" | tee "${sampleID}_fastq_quality_filter.log"
		fastqc "${sampleID}_lengths_and_qual_filtered.fastq"
	# Convert fastq to fasta
		fastq_to_fasta -v -i "${sampleID}_lengths_and_qual_filtered.fastq" -o "${sampleID}_lengths_and_qual_filtered.fasta" | tee "${sampleID}_fastq_to_fasta.log"
	"""
}


// Step 5: Obtain consensus sequences
process consensus_construction {

	publishDir 'Results/5_consensus_sequences', mode: 'copy'
	
	conda "${baseDir}/conda/python_env.yml"

	cpus params.cpus

	input:
	tuple val(sampleID), path(fasta_file)
	path reference_seq_fasta
	val params.cutoff
	val params.threshold

    output:
	path "${sampleID}_consensus.prf"
	path "${sampleID}_discarded.txt"
	path "${sampleID}_frequencies_distribution_graph.png"
	path "${sampleID}_consensus_construction.log"
	tuple val(sampleID), path("${sampleID}_barcodes.json"), emit: barcodes_json_ch
	tuple val("${sampleID}_t${params.threshold}"), path("*_consensus.fna"), env(max_barcode), emit: consensus_reads_ch

	script:
    """
	# Run python script to obtain consensus sequences by barcode
		python ${baseDir}/bin/consensus_construction.py --fasta_file "${fasta_file}" --reference "${reference_seq_fasta}" --output_prefix "${sampleID}" --inferior ${params.cutoff} --threshold_consensus_list ${params.threshold} -c $task.cpus | tee "${sampleID}_consensus_construction.log"
	# Extract the maximum barcode frequency value from the barcodes JSON file
		max_barcode=\$(cat "${sampleID}_barcodes.json" | python -c "import json, sys; data = json.load(sys.stdin); max_val = max(data.values()); print(max_val)")
	"""
}


// Step 6: Map consensus sequences
process consensus_sequences_mapping {

	publishDir 'Results/6_mapped_consensus_sequences', mode: 'copy'
	
	conda "${baseDir}/conda/map_env.yml"
	
	input:
	tuple val(sampleID), path(consensus_fasta), val(max_barcode)
	path reference_seq_fasta
	 
	output:
	path "${sampleID}_align.sam"
	tuple val(sampleID), path("${sampleID}_align_sort.bam"), val(max_barcode), env(consensus_number), emit: bam_file_ch
	
	script:
    """
	# Remove leading 'N' and '-' characters from consensus sequence
		sed -e 's/^N//g' -e 's/-//g' "${consensus_fasta}" > "${sampleID}_filtered_consensus_fasta.fna"
	# Index the reference sequence using BWA
		bwa index "${reference_seq_fasta}"
	# Map consensus sequences to the reference and create a sorted SAM file
		bwa mem "${reference_seq_fasta}" "${sampleID}_filtered_consensus_fasta.fna" | samtools sort -o "${sampleID}_align.sam"
		samtools view -Sb -o "${sampleID}_align.bam" "${sampleID}_align.sam"
		samtools sort "${sampleID}_align.bam" -o "${sampleID}_align_sort.bam"
		samtools index "${sampleID}_align_sort.bam"
	# Count the number of mapped consensus sequences in the SAM file
		consensus_number=\$(awk 'NR > 4 && \$3 != "*" { counter++ } END { print counter }' "${sampleID}_align.sam")
	"""
}


// Step 7: Variant calling
process variant_calling {

	publishDir 'Results/7_vcf_file', mode: 'copy'
	
	conda "${baseDir}/conda/freebayes_env.yml"
	
	input:
	tuple val(sampleID), path(bam_file), val(max_barcode), val(consensus_number)
	path reference_seq_fasta
	
	output:
	tuple val(sampleID), path("${sampleID}.vcf"), val(max_barcode), val(consensus_number),  emit: vcf_ch
	
	script:
    """
	samtools faidx "${reference_seq_fasta}"
	freebayes -p 1 -E 0 -C 1 -F 0 -K --report-monomorphic -dd -f "${reference_seq_fasta}" "${bam_file}" > "${sampleID}.vcf"
	"""
}


// Step 8: VCF analysis
process vcf_analysis {

	publishDir 'Results/8_vcf_analysis', mode: 'copy'
	
	conda "${baseDir}/conda/python_env.yml"
	
	input:
	tuple val(sampleID), path(vcf_file), val(max_barcode), val(consensus_number)
	val params.min_pos
	val params.max_pos
	val params.threshold
	
	output:
	path "*"
	
	script:
    """
	python ${baseDir}/bin/vcf_analyzer.py --input_file "${vcf_file}" --consensus_number "${consensus_number}" --output_prefix "${sampleID}" --min_pos "${params.min_pos}" --max_pos "${params.max_pos}" --max_barcode "${max_barcode}" --cutoff ${params.cutoff} --threshold ${params.threshold}
	"""
}


// Step 9: Offspring search (optional)
process offsprings_finder {

	publishDir 'Results/9_offsprings', mode: 'copy'
	
	conda "${baseDir}/conda/python_env.yml"
	
	input:
	tuple val(sampleID), path(barcodes_json_file)
	
    output:
	path "${sampleID}_differences.png"
	path "${sampleID}_percentage_differences.png"
	path "${sampleID}_offsprings_finder.log"

	script:
    """
	python ${baseDir}/bin/offsprings_finder.py -i "${barcodes_json_file}" --output_prefix "${sampleID}" | tee "${sampleID}_offsprings_finder.log"
	"""
}


// Pipeline
workflow {

	//Step 1: Quality Control
	individual_reads_ch = Channel.fromPath(params.seq_folder_path + '*')
    quality_control(individual_reads_ch)
	
	//Step 2: Joining of paired reads
    reads_pairs_ch = Channel.fromFilePairs("${params.seq_folder_path}*_{1,2}.fastq*")
	reads_merging(reads_pairs_ch)
	
	//Step 3:Graph lengths of the merged reads
	graph_reads_lengths(reads_merging.out.assembled_reads_ch)
	
	// Step 4: Length and quality Filtering | FASTQ to FASTA conversion
	quality_filtering_and_fastq_to_fasta(reads_merging.out.assembled_reads_ch, params.insert_length, params.length_tolerance, params.minimum_quality_score, params.minimum_percent_quality_bases) 
	
	//Step 5: Obtain consensus sequences
	//Value channel "ref_seq_ch" is created to store the path of the reference sequence fasta provided in "params.ref_seq_path"
	ref_seq_ch = Channel.fromPath(params.ref_seq_path).first()
	consensus_construction(quality_filtering_and_fastq_to_fasta.out.fasta_file_ch, ref_seq_ch, params.cutoff, params.threshold)
	
	//Step 6: Map consensus sequences
	consensus_sequences_mapping(consensus_construction.out.consensus_reads_ch, ref_seq_ch)
	
	//Step 7: Variant calling
	variant_calling(consensus_sequences_mapping.out.bam_file_ch, ref_seq_ch)
	
	//Step 8: VCF analysis
	vcf_analysis(variant_calling.out.vcf_ch, params.min_pos, params.max_pos, params.threshold)
	
	//Step 9: Offspring search
	offsprings_finder(consensus_construction.out.barcodes_json_ch)
}


workflow.onComplete {

    println ( workflow.success ? """
        Pipeline summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        WorkDir     : ${workflow.workDir}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}