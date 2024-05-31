# encoding: UTF-8


__doc__ = """
SYNOPSIS
    Python script to find variants in VCF files. 
    It creates an excel file with different data (table with variants, total number of variants, mutation rate...) 
    and 3 graphs showing the distribution of variants in the reference sequence, the distribution of indels and a heatmap with the types of SNPs (if any).

        vcf_analyzer.py -i <file> -c <int> -m <int> -x <int> -o <str> -b <int> -f <int> -t <float> -q <float> [-v,--version] [-h,--help]

DESCRIPTION

 Parameters:
    -i,--input_file	                vcf file [required]
    -c,--consensus_number	        number of consensus sequences used to obtain variants [required]
    -m,--min_pos                    first position of the reference sequence used to quantify mutations [required]
    -x,--max_pos                    last position of the reference sequence used to quantify mutations [required]
    -o,--output_prefix		        output prefix [default='undefined']
    -b,--max_barcode                frequency of the barcode with the maximum frequency [default=None] 
    -f,--cutoff                     lower cut-off of sequences per barcode used previously [default=None]  
    -t,--threshold                  threshold used previously to construct consensus sequences [default=None] 
    -q,--variant_freq_threshold     maximum percentage of consensus sequences with a particular variant (0-100). Variants exceeding this frequency will be excluded from the result [default=95]
                                    variant_freq_threshold=100 will not exclude any variant and variant_freq_threshold=0 will exclude all variants.
    -v,--version	                Shows script version
    -h,--help		                Shows help options
    
    
 Input file:
 The input file must be in VCF format (VCFv4.1) generated by freebayes version 0.9.21.
 
 
 Note: This script was tested using xlsxwriter version 3.1.9. You can install it by running:
    pip install xlsxwriter==3.1.9

AUTHOR

    Javier Martínez del Río    (javier.martinez@cbm.csic.es; javier.mardelrio@gmail.com)
    
"""

__version__ = 'v1.0.0'


## Imports
##----------

import argparse
import re
import sys, os
import pandas as pd
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import math




## Functions
##----------

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input VCF file')
    parser.add_argument('-c', '--consensus_number', type=int, required=True, help='Number of consensus sequences used to obtain variants [required]')
    parser.add_argument('-o', '--output_prefix', type=str, default='undefined', help='Output prefix')
    parser.add_argument('-m', '--min_pos', type=int, required=True, help='First position of the reference sequence used to quantify mutations [required]')
    parser.add_argument('-x', '--max_pos', type=int, required=True, help='Last position of the reference sequence used to quantify mutations [required]')
    parser.add_argument('-b', '--max_barcode', type=int, help='Frequency of the barcode with the maximum frequency')
    parser.add_argument('-f', '--cutoff', type=int, help='Lower cut-off of sequences per barcode used previously ')
    parser.add_argument('-t', '--threshold', type=float, help='Threshold used previously to construct consensus sequences ')
    parser.add_argument('-q', '--variant_freq_threshold', type=float, default=95, help='maximum percentage of consensus sequences with a particular variant. Variants exceeding this frequency will be excluded from the result')
    parser.add_argument("-v","--version", action="version", version=f"%(prog)s {__version__}")
    
    return parser.parse_args()  

def send_error_message(message):
    """Gives an error message"""
    print("\tERROR: {}\n".format(message))
    help()
    sys.exit(1)
         
def convert_cigar_to_string(cigar_string):
    """ 
    Simplifies CIGAR strings to show only letters. Useful for counting types of substitutions in a simpler way.
        
    Parameters:
        - cigar_string (str): The input CIGAR string that represents a sequence alignment. It consists of a series of letters, each with a length, such as "2X2M".

    Returns:
        - simpler_cigar_string (str): A simplified version of the input CIGAR string, where digits and their associated operation letters are expanded. For example, "2X2M" becomes "XXMM".
    """    
    simpler_cigar_string = re.sub(r'(\d+)([MDIX])', lambda match: match.group(2) * int(match.group(1)), cigar_string)
    return simpler_cigar_string

def correct_excessive_matches(row):
    """ 
    Removes redundant letters "M" (matches) that can be generated in CIGARS. For example, a CIGAR of type "2M2X3M" becomes "2X".
    In addition, the convert_cigar_to_string function is used to simplify the CIGAR to "XX". 
    
    Parameters:
        - row (pd.Series): A Pandas Series representing a row of a DataFrame containing information about a genomic variant. 
                      It should have columns 'POS', 'REF', 'ALT', and 'CIGAR'.

    Returns:
        - row (pd.Series): The modified input row after removing excessive matches in the CIGAR string and adjusting the REF and ALT columns accordingly.
    """
    # Simplify the CIGAR string using the convert_cigar_to_string function
    new_cigar = convert_cigar_to_string(row["CIGAR"])
    
    # Count and remove excessive matches at the start of the CIGAR   
    excessive_matches_at_start = 0
    while new_cigar.startswith("M") and not new_cigar[1:].startswith("D") and not new_cigar[1:].startswith("I"):
        new_cigar = new_cigar[1:]
        excessive_matches_at_start += 1
        
    # Count and remove excessive matches at the end of the CIGAR  
    excessive_matches_at_end = 0
    while new_cigar.endswith("M") and not new_cigar[:-1].endswith("D") and not new_cigar[:-1].endswith("I"):
        new_cigar = new_cigar[:-1]
        excessive_matches_at_end += 1
        
    # Adjust the 'POS' column based on matches removed at the start
    row["POS"] = row["POS"] + excessive_matches_at_start
    
    # Adjust the 'REF' and 'ALT' columns based on matches removed at the end
    if excessive_matches_at_end > 0:
        row["REF"] = row["REF"][excessive_matches_at_start:-excessive_matches_at_end]
        row["ALT"] = row["ALT"][excessive_matches_at_start:-excessive_matches_at_end]
    else:
        row["REF"] = row["REF"][excessive_matches_at_start:]
        row["ALT"] = row["ALT"][excessive_matches_at_start:]
        
    # Update the 'CIGAR' column with the modified CIGAR string
    row["CIGAR"] = new_cigar
    return row

def show_one_variant_per_row(dataframe):
    """
    Modifies a DataFrame to represent one variant per row instead of one position per row, as is the case in the original VCF file.

    Parameters:
    - dataframe (pd.DataFrame): The input DataFrame containing genomic variant information. It is assumed to have columns like 'POS', 'REF', 'ALT', 'AO', 'CIGAR', and 'TYPE'.

    Returns:
    - new_dataframe (pd.DataFrame): A new DataFrame where each row corresponds to a specific variant, and information such as 'POS', 'REF', 'ALT', 'AO', 'CIGAR', and 'TYPE' is appropriately duplicated or split based on multiple alternatives.  
    """
    # Initialize an empty list to store rows for the new DataFrame
    new_rows_list = []
    
    # Iterate through each row in the original DataFrame and create a new row for each variant (alternative) from the original row
    for index, row in dataframe.iterrows():
        alternatives = row['ALT'].split(',')
        aos = row['AO'].split(',')
        cigars = row['CIGAR'].split(',')
        types = row['TYPE'].split(',')
        for n, alternative in enumerate(alternatives):
            new_row = [row['POS'], row['REF'], alternatives[n], int(aos[n]), cigars[n], types[n]]
            new_rows_list.append(new_row)
            
    # Create a new DataFrame with the new rows
    new_dataframe = pd.DataFrame(new_rows_list, columns=dataframe.columns)
    return new_dataframe

def create_dataframe_from_vcf(vcf_file, consensus_number, min_pos, max_pos, variant_frequency_threshold):
    """
    Creates a pandas DataFrame from a VCF file, processes and filters the data.

    Parameters:
    - vcf_file (str): Path to the VCF file.
    - consensus_number (int): Number of consensus sequences.
    - min_pos (int): Minimum position to consider.
    - max_pos (int): Maximum position to consider.
    - variant_frequency_threshold (float): maximum percentage of consensus sequences with a particular variant. Variants exceeding this frequency will be excluded from the result

    Returns:
    - df (pd.DataFrame): Processed DataFrame with relevant variant information.
    """
    # Calculate sequence length
    sequence_length = max_pos - min_pos + 1

    # Read VCF file into a DataFrame
    try: 
      df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
      df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown']
    except pd.errors.EmptyDataError:
      print("WARNING: No data found in the VCF file.")
      df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown'])

    # Drop positions with no variants
    df = df[df['ALT'] != "."]

    # Extract relevant information from INFO column
    df['AO'] = df['INFO'].str.extract(r'AO=([\d,]+)').fillna(0) #We extract "AO" from the INFO column (##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">)
    df['CIGAR'] = df['INFO'].str.extract(r'CIGAR=(.*?);') #We extract "CIGAR" from the INFO column (##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">)
    df['TYPE'] = df['INFO'].str.extract(r'TYPE=(.*?)(?:;|$)') #We extract "TYPE" from the INFO column (##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">)

    # Select relevant columns
    df = df[['POS', "REF", "ALT", "AO", "CIGAR", "TYPE"]]

    # Show one variant per row
    df = show_one_variant_per_row(df)

    # Correct excessive matches in CIGAR and update positions
    df = df.apply(correct_excessive_matches, axis=1).sort_values("POS")

    # Filter DataFrame based on min_pos and max_pos
    df = df[(df['POS'] >= min_pos) & (df['POS'] <= max_pos)]

    # Calculate AO rate and variant frequency percentage
    df['AO_rate'] = df['AO'] / (consensus_number * sequence_length) # We add a column with the AO rate of each variant
    df['Variant frequency (%)'] = df['AO'] / consensus_number * 100 # We add a column with the % of consensus sequences with a variant

    # Discard variants with high frequency based on value included in variant_frequency_threshold
    df = df[df['Variant frequency (%)'] < variant_frequency_threshold]

    return df

def create_dict_snp_types(dataframe):
    """
    Creates a dictionary with counts for each type of substitution (SNP variants).

    Parameters:
        - dataframe (pd.DataFrame): Input dataframe created using the create_dataframe_from_vcf function.

    Returns:
        - dict_snp_types (dict): Dictionary containing SNP types as keys and their corresponding counts as values.
    """
    dict_snp_types={}
    def count_snp_types(row):
        snp_type = row["REF"] + row["ALT"]
        try:
            dict_snp_types[snp_type] += row["AO"]
        except:
            dict_snp_types[snp_type] = row["AO"]
    dataframe[dataframe["TYPE"]=="snp"].apply(count_snp_types, axis=1)
    return dict_snp_types

def create_df_snp_types(dict_snp_types):
    """
    Creates a dataframe using the dictionary of snp_types obtained with create_dict_snp_types function.
    In the left part of the table, we show the nucleotides in the reference sequence,
    and in the upper part of the table, the nucleotides for which they are substituted.

    Parameters:
    - dict_snp_types (dict): Dictionary containing counts of SNP types.

    Returns:
    - df_snp_types (pd.DataFrame): DataFrame representing the SNP types table.
    """
    bases = ['G', 'A', 'T', 'C']
    table = {}
    for base in bases:
      table[base] = {}
      for other_base in bases:
          table[base][other_base] = None
    for key, value in dict_snp_types.items():
      base1 = key[1]
      base2 = key[0]
      table[base1][base2] = str(value)
    df_snp_types = pd.DataFrame(table)
    return df_snp_types

def extract_data_from_dataframe(dataframe, output_prefix, consensus_number, min_pos, max_pos, max_barcode, cutoff, threshold):
    """
    Creates an Excel file with different data (table with variants, total number of variants, error rate...), 
    a csv file with the table of variants,
    2 graphs showing the distribution of variants and indel variants,
    and a heatmap graph with the types of SNPs (if any).

    Parameters:
        - dataframe (pd.DataFrame): DataFrame created using create_dataframe_from_vcf function.
        - output_prefix (str): Prefix for output files.
        - consensus_number (int): Number of consensus sequences.
        - min_pos (int): Minimum position analyzed.
        - max_pos (int): Maximum position analyzed.
        - max_barcode (int): Maximum barcode frequency.
        - cutoff (int): Cutoff used.
        - threshold (float): Threshold used.

    Returns:
        - excel_file (str): Path to the Excel file.
        - csv_file (str): Path to the CSV file.
        - variants_graph (str): Path to the variants distribution graph.
        - indels_graph (str): Path to the indels distribution graph.
        - heatmap_snp_types_graph (str): Path to the SNP types heatmap graph.
    """
    #EXTRACT DATA   
    sequence_length = max_pos - min_pos + 1
    dict_snp_types = create_dict_snp_types(dataframe)
    #if all values dict_snp_types are "nan", we assign the value 0 to each substitution type
    if all(math.isnan(value) for value in dict_snp_types.values()): 
        dict_snp_types = {'TC': 0, 'TG': 0, 'AC': 0, 'AG': 0, 'AT': 0, 'GA': 0, 'GC': 0, 'GT': 0, 'CT': 0, 'CA': 0, 'TA': 0, 'CG': 0}
    dataframe = dataframe.dropna()
    total_unique_variants = len(dataframe.index)
    total_variants = dataframe["AO"].sum()
    total_snp = dataframe[dataframe["TYPE"]=="snp"]["AO"].sum()
    total_transitions = sum(dict_snp_types.get(key, 0) for key in ["AG", "GA", "CT", "TC"])
    total_transversions = total_snp - total_transitions
    total_ins = dataframe[dataframe["TYPE"]=="ins"]["AO"].sum()
    total_del = dataframe[dataframe["TYPE"]=="del"]["AO"].sum()
    total_indels = total_ins + total_del
    total_mnp = dataframe[dataframe["TYPE"]=="mnp"]["AO"].sum()
    total_complex = dataframe[dataframe["TYPE"]=="complex"]["AO"].sum()
    try:
        proportion_snp = round((total_snp / total_variants * 100), 2)
        proportion_transitions = round((total_transitions / total_variants * 100), 2)
        proportion_transversions = round((total_transversions / total_variants * 100), 2)
        proportion_indels = round((total_indels / total_variants * 100), 2)
        proportion_ins = round((total_ins / total_variants * 100), 2)
        proportion_del = round((total_del / total_variants * 100), 2) 
        proportion_mnp = round((total_mnp / total_variants * 100), 2)
        proportion_complex = round((total_complex / total_variants * 100), 2)
        mutation_rate = total_variants/(sequence_length*consensus_number)
    except ZeroDivisionError as e:
        print("Error: Cannot divide by zero")
        proportion_snp = 0
        proportion_transitions = 0
        proportion_transversions = 0
        proportion_indels = 0
        proportion_ins = 0
        proportion_del = 0 
        proportion_mnp = 0
        proportion_complex = 0
        mutation_rate = 0
        
    #GROUP DATA IN DICTIONARIES
    basic_data = {
      'Error Rate': [mutation_rate],
      'Total consensus': [consensus_number],
      'Sequence length': [sequence_length],
      'Total Variants': [total_variants],
      'Total Unique Variants': [total_unique_variants],
      'Max barcode frequency': [max_barcode],
      'Cutoff used': [cutoff],
      'Threshold used': [threshold],
      'Min position analysed': [min_pos],
      'Max position analysed': [max_pos],     
    }
    total_data = {
      'Total Variants': [total_variants],
      'Total SNP': [total_snp],
      'Total Transitions': [total_transitions],
      'Total Transversions': [total_transversions],
      'Total INS': [total_ins],
      'Total DEL': [total_del],
      'Total Indels': [total_indels],
      'Total MNP': [total_mnp],
      'Total Complex': [total_complex]
    }
    proportion_data = {
      'Proportion SNP': [proportion_snp],
      'Proportion SNP Transitions': [proportion_transitions],
      'Proportion SNP Transversions': [proportion_transversions],
      'Proportion Indels': [proportion_indels],
      'Proportion INS': [proportion_ins],
      'Proportion DEL': [proportion_del],
      'Proportion MNP': [proportion_mnp],
      'Proportion complex': [proportion_complex]
    }

    #CREATE DATAFRAMES
    df_grouped_by_position = dataframe.groupby('POS')[['AO', 'AO_rate', 'Variant frequency (%)']].sum().reset_index().sort_values('POS')
    df_indels = dataframe[(dataframe["TYPE"] == "ins") | (dataframe["TYPE"] == "del")]
    df_grouped_indels = df_indels.groupby('POS')[['AO', 'AO_rate', 'Variant frequency (%)']].sum().reset_index().sort_values('POS') # Group indels' AO values by position and sum them up
    df_snps = dataframe[(dataframe["TYPE"] == "snp")]
    df_grouped_snp = df_snps.groupby('POS')[['AO', 'AO_rate', 'Variant frequency (%)']].sum().reset_index().sort_values('POS')
    df_snp_types = create_df_snp_types(dict_snp_types)  
    df_basic_data = pd.DataFrame(basic_data)
    df_total_data = pd.DataFrame(total_data)
    df_proportion_data = pd.DataFrame(proportion_data)


    #CREATE EXCEL FILE FROM DATAFRAMES
    excel_file = output_prefix + '.xlsx'
    with pd.ExcelWriter(excel_file, engine='xlsxwriter') as writer:
        dataframe.to_excel(writer, sheet_name='Raw data', index=False)
        df_grouped_by_position.to_excel(writer, sheet_name='Variants per position', index=False)
        df_indels.to_excel(writer, sheet_name='Indels', index=False)
        df_grouped_indels.to_excel(writer, sheet_name='Indels per position', index=False)
        df_snps.to_excel(writer, sheet_name='SNPs', index=False)
        df_grouped_snp.to_excel(writer, sheet_name='SNPs per position', index=False)
        df_snp_types.to_excel(writer, sheet_name='SNP types', index=True)
        df_basic_data.to_excel(writer, sheet_name='Error rate', index=False)
        df_total_data.to_excel(writer, sheet_name='Total data', index=False)
        df_proportion_data.to_excel(writer, sheet_name='Proportion data', index=False)


    #CREATE CSV FILE FROM DATAFRAME WITH RAW DATA
    csv_file = output_prefix + ".csv"
    dataframe.to_csv(csv_file, index=False)


    #CREATE GRAPH WITH THE DISTRIBUTION OF VARIANTS
    plt.clf()
    try:
      grouped_df = dataframe.groupby('POS')['AO'].sum().reset_index().sort_values('POS') # Group AO values by position and sum them up
    except KeyError:
        print("WARNING: There is a problem with grouped_df Dataframe. The variants distribution graph will be empty.") 
    plt.bar(grouped_df['POS'], grouped_df['AO'], color='tab:blue', label='Variants')
    plt.xlabel('Position')
    plt.ylabel('Number of variants')
    plt.title('Variants distribution')
    plt.suptitle(output_prefix)
    plt.legend(loc='upper right')
    try:
        plt.ylim([0, grouped_df['AO'].max()*1.1]) # Set y-axis limits
    except Exception:
        pass
    plt.xlim(min_pos -1, max_pos + 1) # Set x-axis limits
    variants_graph = output_prefix + "_variants_distribution.png"
    plt.savefig(variants_graph)
    plt.close()


    #CREATE GRAPH WITH THE DISTRIBUTION OF INDEL VARIANTS
    plt.clf()
    try:
      plt.bar(df_grouped_indels['POS'], df_grouped_indels['AO'], color='tab:red', label='Indel variants')
    except KeyError:
        print("WARNING: There is a problem with df_grouped_indels Dataframe. The indels graph will be empty.")
    plt.xlabel('Position')
    plt.ylabel('Number of variants')
    plt.title('Indel variants distribution')
    plt.suptitle(output_prefix)
    plt.legend(loc='upper right')
    try:
        plt.ylim([0, df_grouped_indels['AO'].max()*1.1]) # Set y-axis limits
    except Exception:
        pass
    plt.xlim(min_pos -1, max_pos + 1) # Set x-axis limits
    # Save graph
    indels_graph = output_prefix + "_indels_distribution.png"
    plt.savefig(indels_graph)
    plt.close()


    #CREATE HEATMAP IMAGE OF THE SNP TYPES (shown as percentages)
    plt.clf()
    total_sum = sum(dict_snp_types.values())
    if total_sum != 0:
        dict_snp_types_percentage = {key: (value / total_sum) * 100 for key, value in dict_snp_types.items()}
    else:
        dict_snp_types_percentage = dict_snp_types
    df_table = create_df_snp_types(dict_snp_types_percentage)
    plt.figure(figsize=(6, 4))
    sns.heatmap(df_table.astype(float), cmap='coolwarm', annot=True, fmt=".1f", cbar_kws={'label': '%'})
    plt.gca().xaxis.set_label_position('top')
    plt.xlabel('Nucleotide substitutions')
    plt.ylabel('Reference nucleotides')
    plt.suptitle(output_prefix)
    plt.subplots_adjust(top=0.80)
    plt.gca().tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False) #Adjust position of ticks and x-axis labels
    plt.text(0.5, -0.1, f"Total number of substitutions: {total_sum}", transform=plt.gca().transAxes, ha='center') #Add total number of substitutions as text in the image
    heatmap_snp_types_graph = output_prefix + 'heatmap_snp_types.png'
    plt.savefig(heatmap_snp_types_graph)
    plt.close()
    
    return excel_file, csv_file, variants_graph, indels_graph, heatmap_snp_types_graph


## Main program
##-------------

def main():


    ## Step 1: Parameters catching
    ##-------------------- 

    print("\tStep 1 => Parameters catching")

    args = parse_arguments()
    vcf_file = args.input_file
    consensus_number = args.consensus_number
    output_prefix = args.output_prefix
    min_pos = args.min_pos
    max_pos = args.max_pos
    max_barcode = args.max_barcode
    cutoff = args.cutoff
    threshold = args.threshold
    variant_frequency_threshold = args.variant_freq_threshold 
   
    if not os.path.isfile(vcf_file):
        send_error_message(f'Input file "{vcf_file}" of {__file__} does not exist. Please provide a valid file path.')
        
    print("\tStep 1 done\n")



    ## Step 2: Creating dataframe from vcf file
    ##-------------------- 

    print("\tStep 2 => Creating dataframe from vcf file")
    
    try:
        dataframe_from_vcf = create_dataframe_from_vcf(vcf_file, consensus_number, min_pos, max_pos, variant_frequency_threshold)
    except Exception as e:
        send_error_message(f'VCF file could not be analysed: Problem when trying to convert VCF file to pandas dataframe:\n{e}') 
    
    print("\tStep 2 done\n")



    ## Step 3: Extracting data from dataframe | Creating excel and csv files | Drawing graphs
    ##-------------------- 
    
    print("\tStep 3 => Extracting data from dataframe | Creating excel and csv files | Drawing graphs")
    try:
        excel_file, csv_file, variants_graph, indels_graph, heatmap_snp_types_graph = extract_data_from_dataframe(dataframe_from_vcf, output_prefix, consensus_number, min_pos, max_pos, max_barcode, cutoff, threshold)
    except Exception as e:
        send_error_message(f'There was a problem with the extract_data_from_dataframe function:\n{e}') 
    
    print("\tStep 3 done\n")
    
    print("\tJOB DONE!")

if __name__ == "__main__":
    main()