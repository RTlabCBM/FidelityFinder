# encoding: UTF-8


# Imports
import re
import os
import sys
import pandas as pd
import xlsxwriter
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import math


__doc__ = """
SYNOPSIS
    Python script to find variants in VCF files. Creates an excel file with different data (table with variants, total number of variants, mutation rate...) and 3 graphs showing the distribution of variants in the reference sequence, the distribution of indels and a heatmap with the types of SNPs (if any).

        vcf_analyzer.py -i <file> -c <int> -m <int> -x <int> -o <str> [-h,--help]

DESCRIPTION

 Parameters:
    -i,--input_vcf	            vcf file [required]
    -c,--consensus_number	    number of consensus sequences used to obtain variants [required]
    -o,--output_prefix		    output prefix [default='undefined']
    -m,--min_pos                first position of the reference sequence used to quantify mutations [required]
    -x,--max_pos                last position of the reference sequence used to quantify mutations [required]
    -h,--help		            Shows help options


AUTHOR

    Javier Martínez del Río    <javier.martinez@cbm.csic.es>
    
"""


## Functions
##----------

def help():
    print(globals()['__doc__'])
    sys.exit(1)

def show_error_message(message):
    """Gives an error message"""
    print("\tERROR: {}\n".format(message))
    help()
    sys.exit(1)
    
def show_error_messageB(param, message):
    """Gives an error message if argument is not valid"""
    print("\tERROR: '{} {}'\n".format(param, message))
    help()
    sys.exit(1)      

def convert_cigar_to_string(cadena):
    # Simplifies CIGAR strings to show only letters. For example, converts CIGAR "1X" to "X" or CIGAR "2X2M" to "XXMM".  Useful for counting type of substitutions in a simpler way.
    resultado = re.sub(r'(\d+)([MDIX])', lambda match: match.group(2) * int(match.group(1)), cadena)
    return resultado

def correct_excessive_matches(row):
    # Removes redundant letters "M" (matches) that can be generated in CIGARS. For example, a CIGAR of type "2M2X3M" becomes "2X". In addition, the convert_cigar_to_string function is used to simplify the CIGAR to "XX". When redudant letters "M" are found, the data in the REF and ALT columns are also corrected accordingly.
    new_cigar = convert_cigar_to_string(row["CIGAR"])
    excessive_matches_at_start = 0
    excessive_matches_at_end = 0
    while new_cigar.startswith("M") and not new_cigar[1:].startswith("D") and not new_cigar[1:].startswith("I"):
        new_cigar = new_cigar[1:]
        excessive_matches_at_start += 1
    while new_cigar.endswith("M") and not new_cigar[:-1].endswith("D") and not new_cigar[:-1].endswith("I"):
        new_cigar = new_cigar[:-1]
        excessive_matches_at_end += 1
    row["POS"] = row["POS"] + excessive_matches_at_start
    if excessive_matches_at_end > 0:
        row["REF"] = row["REF"][excessive_matches_at_start:-excessive_matches_at_end]
        row["ALT"] = row["ALT"][excessive_matches_at_start:-excessive_matches_at_end]
    else:
        row["REF"] = row["REF"][excessive_matches_at_start:]
        row["ALT"] = row["ALT"][excessive_matches_at_start:]
    row["CIGAR"] = new_cigar
    return row

def show_one_variant_per_row(dataframe):
    #Modifies dataframe to show one variant per row, instead of one position per row as is the case in the original VCF file
    new_rows = []
    for index, row in dataframe.iterrows():
        alternatives = row['ALT'].split(',')
        aos = row['AO'].split(',')
        cigars = row['CIGAR'].split(',')
        types = row['TYPE'].split(',')
        for n, alternative in enumerate(alternatives):
            new_row = [row['POS'], row['REF'], alternatives[n], int(aos[n]), cigars[n], types[n]]
            new_rows.append(new_row)
    new_dataframe = pd.DataFrame(new_rows, columns=dataframe.columns)
    return new_dataframe

def create_dataframe_from_vcf(vcf_file, consensus_number, min_pos, max_pos):
  sequence_length = max_pos - min_pos + 1
  # Create a pandas dataframe from the input vcf file. We also need the number of consensus and the length of the sequence to calculate the rate of the variants.
  df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
  df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown']
  df = df[df['ALT'] != "."] #We drop positions with no variants
  if len(df) < 1:
    print("VCF file has no variants")
  df['AO'] = df['INFO'].str.extract(r'AO=([\d,]+)').fillna(0) #We extract "AO" from the INFO column (##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">)
  df['CIGAR'] = df['INFO'].str.extract(r'CIGAR=(.*?);') #We extract "CIGAR" from the INFO column (##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">)
  df['TYPE'] = df['INFO'].str.extract(r'TYPE=(.*?)(?:;|$)') #We extract "TYPE" from the INFO column (##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">)
  df = df[['POS', "REF", "ALT", "AO", "CIGAR", "TYPE"]]
  df = show_one_variant_per_row(df)
  df = df.apply(correct_excessive_matches, axis=1).sort_values("POS")
  df = df[(df['POS'] >= min_pos) & (df['POS'] <= max_pos)] #We filter the df according to the input parameters "min_pos" and "max_pos" 
  df['AO_rate'] = df['AO'] / (consensus_number * sequence_length) # We add a column with the AO rate of each variant
  return df

def create_dict_snp_types(dataframe):
  # Creates a dictionary with each type of substitution (SNP variants). Input dataframe: dataframe created using create_dataframe_from_vcf function.
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
  # Creates a dataframe using the dictionary of snp_types obtained with create_dict_snp_types function
  # In the left part of the table we show the nucleotides in the reference sequence and in the upper part of the table the nucleotides for which they are substituted.
  bases = ['G', 'A', 'T', 'C']
  table = {}
  for base in bases:
      table[base] = {}
      for otra_base in bases:
          table[base][otra_base] = None
  for key, value in dict_snp_types.items():
      base1 = key[1]
      base2 = key[0]
      table[base1][base2] = str(value)
  df_snp_types = pd.DataFrame(table)
  return df_snp_types

#def find_numeric_columns(df):
#    numeric_columns  = df.select_dtypes(include=['int64', 'float64']).columns
#    return numeric_columns 

def extract_data_from_dataframe(dataframe, output_file_name, consensus_number, min_pos, max_pos):
  # Creates an excel file with different data (table with variants, total number of variants, mutation rate...) and 3 graphs showing the distribution of variants in the reference sequence, the distribution of indels and a heatmap with the types of SNPs (if any).
  # Input dataframe: dataframe created using create_dataframe_from_vcf function.

  #Extract data
  sequence_length = max_pos - min_pos + 1
  dict_snp_types = create_dict_snp_types(dataframe)
  if all(math.isnan(value) for value in dict_snp_types.values()): #if all values dict_snp_types are "nan", we assign the value 0 to each substitution type
    dict_snp_types = {'TC': 0, 'TG': 0, 'AC': 0, 'AG': 0, 'AT': 0, 'GA': 0, 'GC': 0, 'GT': 0, 'CT': 0, 'CA': 0, 'TA': 0, 'CG': 0}
  #total_nan_values = dataframe.isna().sum().sum()
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
        
  #Group data in dictionaries
  basic_data = {
      'Total consensus': [consensus_number],
      'Sequence length': [sequence_length],
      'Total Variants': [total_variants],
      'Total Unique Variants': [total_unique_variants],
      'Mutation Rate': [mutation_rate]
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

  #Create dataframes
  df_grouped_by_position = dataframe.groupby('POS')[['AO', 'AO_rate']].sum().reset_index().sort_values('POS')
  df_indels = dataframe[(dataframe["TYPE"] == "ins") | (dataframe["TYPE"] == "del")]
  df_grouped_indels = df_indels.groupby('POS')[['AO', 'AO_rate']].sum().reset_index().sort_values('POS') # Group indels' AO values by position and sum them up
  df_snps = dataframe[(dataframe["TYPE"] == "snp")]
  df_grouped_snp = df_snps.groupby('POS')[['AO', 'AO_rate']].sum().reset_index().sort_values('POS')
  df_snp_types = create_df_snp_types(dict_snp_types)  
  df_basic_data = pd.DataFrame(basic_data)
  df_total_data = pd.DataFrame(total_data)
  df_proportion_data = pd.DataFrame(proportion_data)

  
  # Create excel file from dataframes
  excel_file = output_file_name + '.xlsx'
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
  

  # Save dataframe with raw data as csv file
  csv_file = output_file_name + ".csv"
  dataframe.to_csv(csv_file, index=False)


  #Create heatmap image of the SNP types (shown as percentages)
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
  plt.suptitle(output_file_name)
  plt.subplots_adjust(top=0.80)
  plt.gca().tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False) #Adjust position of ticks and x-axis labels
  plt.text(0.5, -0.1, f"Total number of substitutions: {total_sum}", transform=plt.gca().transAxes, ha='center') #Add total number of substitutions as text in the image
  name_heatmap_snp_types = output_file_name + 'heatmap_snp_types.png'
  plt.savefig(name_heatmap_snp_types)
  plt.close()

  #Create graph with the distribution of variants in the reference sequence
  plt.clf()
  grouped_df = dataframe.groupby('POS')['AO'].sum().reset_index().sort_values('POS') # Group AO values by position and sum them up
  colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
  plt.bar(grouped_df['POS'], grouped_df['AO'], color=colors[0])
  plt.xlabel('Position')
  plt.ylabel('Number of variants')
  plt.title('Variants distribution')
  plt.suptitle(output_file_name)
  mean = grouped_df['AO'].mean()
  plt.axhline(y=mean, color=colors[1], linestyle='--', label=f'Mean ({mean:.2f})') # Add a horizontal line to mark the mean of variants in positions with variants
  median = grouped_df['AO'].median()
  plt.axhline(y=median, color=colors[2], linestyle='--', label=f'Median ({median:.2f})') # Add a horizontal line to mark the median of variants in positions with variants
  plt.legend(['Mean ({:.2f})'.format(mean), 'Median ({:.2f})'.format(median), 'Variants'], loc='upper right')
  try:
    plt.ylim([0, grouped_df['AO'].max()*1.1]) # Set y-axis limits
  except Exception:
    pass
  plt.xlim(min_pos -1, max_pos + 1) # Set x-axis limits
  name_variants_graph = output_file_name + "_variants_distribution.png"
  plt.savefig(name_variants_graph)
  plt.close()

  #Create graph with the distribution of indel variants in the reference sequence
  plt.clf()
  colors = ['red', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
  plt.bar(df_grouped_indels['POS'], df_grouped_indels['AO'], color=colors[0])
  plt.xlabel('Position')
  plt.ylabel('Number of variants')
  plt.title('Indel variants distribution')
  plt.suptitle(output_file_name)
  mean = df_grouped_indels['AO'].mean()
  plt.axhline(y=mean, color=colors[1], linestyle='--', label=f'Mean ({mean:.2f})') # Add a horizontal line to mark the mean of variants in positions with variants
  median = df_grouped_indels['AO'].median()
  plt.axhline(y=median, color=colors[2], linestyle='--', label=f'Median ({median:.2f})') # Add a horizontal line to mark the median of variants in positions with variants
  plt.legend(['Mean ({:.2f})'.format(mean), 'Median ({:.2f})'.format(median), 'Indels'], loc='upper right')
  try:
    plt.ylim([0, df_grouped_indels['AO'].max()*1.1]) # Set y-axis limits
  except Exception:
    pass
  plt.xlim(min_pos -1, max_pos + 1) # Set x-axis limits
  # Descargar la gráfica
  indels_graph = output_file_name + "_indels_distribution.png"
  plt.savefig(indels_graph)
  plt.close()


## Parameters catching
##--------------------

# defaults
output_prefix = "undefined"

print("\tStep 1 => Parameters catching")
script_name = sys.argv.pop(0)
if len(sys.argv) == 0:
  show_error_message("You have to provide some parameters")
else:
  while len(sys.argv) > 0:
        param = sys.argv.pop(0)
        if param == '-i' or param == '--input_vcf':
            vcf_file = sys.argv.pop(0)
            if not os.path.exists(vcf_file):
                show_error_messageB(vcf_file, "file does not exist !!!")
        elif param == '-c' or param == '--consensus_number':
            consensus_number = int(sys.argv.pop(0))
        elif param == '-o' or param == '--output_prefix':
            output_prefix = sys.argv.pop(0)
        elif param == '-m' or param == '--min_pos':
            min_pos = int(sys.argv.pop(0))
        elif param == '-x' or param == '--max_pos':
            max_pos = int(sys.argv.pop(0))            
        elif param == '-h' or param == '--help':
            help()
            sys.exit(1)
        else:
            if param == script_name:
                continue
            else:
                show_error_messageB(script_name, param, "parameter not recognized !!!")
print("\tStep 1 done\n")




## Main program
##-------------

def main():

    print("\tStep 2 => Creating dataframe")
    try:
        dataframe_from_vcf = create_dataframe_from_vcf(vcf_file, consensus_number, min_pos, max_pos)
    except Exception as e:
        show_error_message(f'VCF file could not be analysed: Problem when trying to convert VCF file to pandas dataframe:\n{e}') 
    print("\tStep 2 done\n")

    print("\tStep 3 => Extracting data from dataframe and drawing graphs")
    try:
        extract_data_from_dataframe(dataframe_from_vcf, output_prefix, consensus_number, min_pos, max_pos)
    except Exception as e:
        show_error_message(f'There was a problem with the extract_data_from_dataframe function:\n{e}') 
    print("\tStep 3 done\n")
    
    print("\tJOB DONE!")

if __name__ == "__main__":
    main()