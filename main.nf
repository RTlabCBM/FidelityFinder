#!/usr/bin/env nextflow


/* 

	############################################################
	###################### FidelityFinder ######################
	############################################################
	
	DESCRIPTION
	FidelityFinder is a bioinformatics analysis pipeline to determine the fidelity transcriptases and the fidelity of DNA synthesis by 
	reverse transcriptases (RTs) from sequences obtained by Next Generation Sequencing (NGS).
	
	AUTHOR
	Javier Martínez del Río (javier.martinez@cbm.csic.es)
	
	LICENSE
	Attribution-NonCommercial-ShareAlike 4.0 International
	(CC BY-NC-SA 4.0) (https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en)

	FidelityFinder includes several third-party packages provided under other open source licenses, please check them for additional details.
	
*/




nextflow.enable.dsl=2

// Step 1: Quality Control using FASTQC (optional)
process fastqc {

	publishDir 'Results/1_fastqc_results', mode: 'copy'

	conda "${baseDir}/ff_env.yml"
	
    input:
    path fastq
	
	output:
	path "*"
	
    script:
    """
    fastqc ${fastq}
    """
}

// Step 2: Merging of reads with PEAR
process pear {

    publishDir 'Results/2_merged_sequences', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
    
	input:
	tuple val(sampleId), file(reads)
    val params.insert_length
	
	output:   
    path "${sampleId}.unassembled.forward.fastq"
    path "${sampleId}.unassembled.reverse.fastq"
    path "${sampleId}.discarded.fastq"
	path "${sampleId}.assembled.fastq"
    path "${sampleId}_pear.log"
	tuple val(sampleId), path("${sampleId}.assembled_filtered.fastq"), emit: assembled_seq_ch
	tuple val(sampleId), path("${sampleId}_lengths.txt"), emit: len_ch
	

	script:
    """
    pear -f ${reads[0]} -r ${reads[1]} -j ${task.cpus} -o "${sampleId}" | tee "${sampleId}_pear.log" 
	awk 'NR % 4 == 2 {dict_lengths[length(\$0)]++} END {for (l in dict_lengths) {print l, dict_lengths[l]}}' "${sampleId}.assembled.fastq" | sort -n > "${sampleId}_lengths.txt" 
	awk 'BEGIN {FS = "\\t" ; OFS = "\\n"} {header = \$0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= ("${params.insert_length}"-20) && length(seq) <= ("${params.insert_length}"+20)) {print header, seq, qheader, qseq}}' "${sampleId}.assembled.fastq" > "${sampleId}.assembled_filtered.fastq"  
	"""
}


// Step 3: Graph lengths of the merged reads
process graph_lengths {

	publishDir 'Results/3_len_graphs', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val len_file
	
    output:
	path "*"
	
	script:
    """
	python ${baseDir}/bin/sizes_graphs.py --input_file "${len_file[1]}" --output_prefix "${len_file[0]}"
	"""
}


// Step 4: Obtain sequences in fasta format
process fastq_to_fasta {

	publishDir 'Results/4_fasta sequences', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val assembled_seq_ch
	
    output:
	path "${assembled_seq_ch[0]}_qual_filtered.fastq"
	path "${assembled_seq_ch[0]}_fastq_quality_filter.log"
	tuple val("${assembled_seq_ch[0]}"), path("${assembled_seq_ch[0]}_qual_filtered.fasta"), emit: fasta_ch
	path "${assembled_seq_ch[0]}_fastq_to_fasta.log"
	
	script:
    """	
	fastq_quality_filter -v -q 1 -p 60 -i "${assembled_seq_ch[1]}" -o "${assembled_seq_ch[0]}_qual_filtered.fastq" | tee "${assembled_seq_ch[0]}_fastq_quality_filter.log"
	fastq_to_fasta -v -i "${assembled_seq_ch[0]}_qual_filtered.fastq" -o "${assembled_seq_ch[0]}_qual_filtered.fasta" | tee "${assembled_seq_ch[0]}_fastq_to_fasta.log"
	"""
}

// Step 5: Obtain consensus sequences
process consensus_construction {

	publishDir 'Results/5_consensus', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val fasta_ch
	path ref_seq_path
	val params.fw_primer
	val params.rv_primer
	val params.cutoff
	val params.bc_size
	val params.threshold

    output:
	path "${fasta_ch[0]}_consensus.fna"
	path "${fasta_ch[0]}_consensus.png"
	path "${fasta_ch[0]}_consensus.prf"
	path "${fasta_ch[0]}_consensus.xls"
	path "${fasta_ch[0]}_cutoff_consensus.png"
	path "${fasta_ch[0]}_discarded.txt"
	tuple val("${fasta_ch[0]}"), path("${fasta_ch[0]}_barcodes.json"), emit: barcodes_ch
	tuple val("${fasta_ch[0]}"), path("${fasta_ch[0]}_filtered_consensus_without_gaps.fna"), emit: consensus_ch
	path "${fasta_ch[0]}_consensus_construction.log"
	
	script:
    """
	python ${baseDir}/bin/consensus_construction_barcode_finder.py -f "${fasta_ch[1]}" --reference "${ref_seq_path}" --output_prefix "${fasta_ch[0]}" --fw_primer ${params.fw_primer}  --rv_primer ${params.rv_primer} --inferior ${params.cutoff} --bc_size ${params.bc_size} --threshold_consensus ${params.threshold} | tee "${fasta_ch[0]}_consensus_construction.log"
	sed 's/^N//g' "${fasta_ch[0]}_consensus.fna" > "${fasta_ch[0]}_filtered_consensus.fna"
	sed 's/-//g' "${fasta_ch[0]}_filtered_consensus.fna" > "${fasta_ch[0]}_filtered_consensus_without_gaps.fna"
	"""
}


// Step 6: Map consensus sequences
process consensus_sequences_mapping {

	publishDir 'Results/6_mapped_consensus_sequences', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val consensus_ch
	path ref_seq_path
	 
	output:
	path "${consensus_ch[0]}_align.sam"
	tuple val("${consensus_ch[0]}"), path("${consensus_ch[0]}_align_sort.bam"), env(consensus_number), emit: bam_file_ch
	env avg_seq_len, emit: avg_seq_len_ch //this is no longer needed
	
	script:
    """
	bwa index "${ref_seq_path}"
	bwa mem "${ref_seq_path}" "${consensus_ch[1]}" | samtools sort -o "${consensus_ch[0]}_align.sam"
	samtools view -Sb -o "${consensus_ch[0]}_align.bam" "${consensus_ch[0]}_align.sam"
	samtools sort "${consensus_ch[0]}_align.bam" -o "${consensus_ch[0]}_align_sort.bam"
	samtools index "${consensus_ch[0]}_align.bam"
	consensus_number=\$(awk 'NR > 4 && \$3 != "*" { counter++ } END { print counter }' "${consensus_ch[0]}_align.sam")
	avg_seq_len=\$(awk 'NR > 4 && \$3 != "*" { sum += length(\$10); count++; } END { average = sum / count; print average }' "${consensus_ch[0]}_align.sam") 
	"""
}


// Step 7: Variant calling
process variant_calling {

	publishDir 'Results/7_vcf_file', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val bam_file_ch
	path ref_seq_path
	
	output:
	tuple val("${bam_file_ch[0]}"), path("${bam_file_ch[0]}.vcf"), val("${bam_file_ch[2]}"),  emit: vcf_ch
	
	script:
    """
	samtools faidx "${ref_seq_path}"
	freebayes -p 1 -E 0 -C 1 -F 0 -K --report-monomorphic -f "${ref_seq_path}"  "${bam_file_ch[1]}" > "${bam_file_ch[0]}.vcf"
	"""
}


// Step 8: VCF analysis
process vcf_analysis {

	publishDir 'Results/8_vcf_analysis', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val vcf_ch
	val params.min_pos
	val params.max_pos
	
	output:
	path "*"
	
	script:
    """
	python ${baseDir}/bin/vcf_analyzer.py --input_vcf "${vcf_ch[1]}" --consensus_number "${vcf_ch[2]}" --output_prefix "${vcf_ch[0]}" --min_pos "${params.min_pos}" --max_pos "${params.max_pos}"
	"""
}


// Step 9: Offspring search (optional)
process offsprings_finder {

	publishDir 'Results/9_offsprings', mode: 'copy'
	
	conda "${baseDir}/ff_env.yml"
	
	input:
	val barcodes_ch
	
    output:
	path "*"
	
	script:
    """
	python ${baseDir}/bin/offsprings_finder.py -f "${barcodes_ch[1]}" --output_prefix "${barcodes_ch[0]}" | tee "${barcodes_ch[0]}_offsprings_finder.log"
	"""
}




// Pipeline
workflow {
	fastq_ch = Channel.fromPath(params.seq_folder_path + '*')
    fastqc(fastq_ch)
    pear_ch = Channel.fromFilePairs("${params.seq_folder_path}*_{1,2}.fastq*")
	pear(pear_ch, params.insert_length)
	graph_lengths(pear.out.len_ch)
	//fastq_to_fasta(pear.out.assembled_seq_ch) 
	//ref_seq_ch = Channel.fromPath(params.ref_seq_path).first()
	//consensus_construction(fastq_to_fasta.out.fasta_ch, ref_seq_ch, params.fw_primer, params.rv_primer, params.cutoff, params.bc_size, params.threshold)
	//consensus_sequences_mapping(consensus_construction.out.consensus_ch, ref_seq_ch)
	//variant_calling(consensus_sequences_mapping.out.bam_file_ch, ref_seq_ch)
	//vcf_analysis(variant_calling.out.vcf_ch, params.min_pos, params.max_pos)
	//offsprings_finder(consensus_construction.out.barcodes_ch)
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
