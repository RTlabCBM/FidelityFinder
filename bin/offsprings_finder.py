# encoding:UTF-8


# Imports
from Bio import SeqIO
import json
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys, os, math
from collections import defaultdict


########################################################################################

__doc__ = """
SYNOPSIS

 Python script to find possible offspring barcodes (barcodes generated due to errors in PCR reactions/sequencing reaction).
 This script finds barcodes that have only one or two different nucleotides with respect to other barcodes with a frequency equal or higher.

    offsprings_finder.py -f <file> [-o,--output] [-h,--help] [-v,--version]

DESCRIPTION

 Parameters:
    -f,--file		        json file [required]
    -o,--output_prefix		output prefix [default='undefined']
    -v,--version	        Shows script version
    -h,--help		        Shows help options


AUTHORS

    Javier Martinez del Río (javier.martinez@cbm.csic.es)
                    
"""
__author__ = "Javier Martínez del Río"
__version__ = 'v1'

## Functions
##----------

def help():
    print(globals()['__doc__'])
    sys.exit(1)

def paramerror(param, message):
    """Gives an error message if argument is not valid"""
    print("\tERROR: '{} {}'\n".format(param, message))
    help()
    sys.exit(1)

def paramerrorB(message):
    """Gives an error message if argument is not valid"""
    print("\tERROR: {}\n".format(message))
    help()
    sys.exit(1)

def version():
    """Print the version of this script"""
    print("\n{}".format(script_name))
    print("Version: {}".format(__version__))
    print("Authors: {}\n".format(__author__))
    sys.exit(1)


## Parameters catching
##--------------------


# defaults
out_prefix    = "undefined"


print("\tStep 1 => Parameters catching")
script_name = sys.argv.pop(0)
if len(sys.argv) == 0:
  paramerrorB("You have to provide some parameters !!!")
  version()
else:
  while len(sys.argv) > 0:
        param = sys.argv.pop(0)
        if param == '-f' or param == '--barcodes_file':
            barcodes_file = sys.argv.pop(0)
            if not os.path.exists(barcodes_file):
                paramerror(barcodes_file, "file does not exist !!!")
        elif param == '-o' or param == '--output_prefix':
            out_prefix = sys.argv.pop(0)               
        elif param == '-h' or param == '--help':
            help()
            sys.exit(1)
        elif param == '-v' or param == '--version':
            version()
            sys.exit(1)
        else:
            if param == script_name:
                continue
            else:
                paramerror(param, "parameter not recognized !!!")
print("\tStep 1 done\n")
print("\n")





## Main program
##-------------

print("\tStep 2 => Open & read input file")
with open(barcodes_file, 'r') as json_file:
    barcodes_dic = json.load(json_file)
print("\tStep 2 done\n")


print("\tStep 3 => Finding differences between barcodes")
dic_1dif={} #we store as keys barcodes that differ in 1 nt from other barcodes, the value of each key is a list with the barcodes with 1 difference
dic_2dif={} #we store as keys barcodes that differ in 2 nt from other barcodes, the value of each key is a list with the barcodes with 2 differences
barcodes_list = barcodes_dic.keys()
for seq1 in barcodes_list:
  for seq2 in barcodes_list:
      differences_count = 0
      for nt in range(len(seq1)):
        if seq1[nt] != seq2[nt]:
          differences_count +=1
          if differences_count == 3:
            break
      if differences_count == 1:         
        try:
          dic_1dif[seq1].append(seq2)
        except KeyError:
          dic_1dif[seq1]=[seq2]
      elif differences_count == 2:         
        try:
          dic_2dif[seq1].append(seq2)
        except KeyError:
          dic_2dif[seq1]=[seq2]
print("\t\tBarcodes with 1 difference from other barcodes: ", len(dic_1dif))
print("\t\tBarcodes with 2 differences from other barcodes: ", len(dic_2dif))
print("\tStep 3 done\n")




print("\tStep 4 => Finding frequencies of barcodes with 1 or 2 differences from other barcodes with equal or higher frequency")

#We create the dict "dic_freq_1dif", the keys are frequencies of barcodes and the value of each key is the number of barcodes of that frequency that differ in only one nt from other barcodes of equal or higher frequency
dic_freq_1dif={}
for barcode,similar_barcodes in dic_1dif.items():
  barcode_freq = barcodes_dic[barcode]
  list_freq_values = [barcodes_dic[similar_bc] for similar_bc in similar_barcodes] 
  if barcode_freq <= max(list_freq_values): 
    try:
      dic_freq_1dif[barcode_freq] +=1
    except:
      dic_freq_1dif[barcode_freq] =1
dic_freq_1dif_list = sorted(list(dic_freq_1dif.items())) #We order dic_freq_1dif by frequencies:
dic_freq_1dif = dict(dic_freq_1dif_list)
print("\t\tFrequencies of barcodes with 1 difference from other barcodes with equal or higher frequency: ", dic_freq_1dif)

#We create the dict "dic_freq_2dif", the keys are frequencies of barcodes and the value of each key is the number of barcodes of that frequency that differ in two nt from other barcodes of equal or higher frequency
dic_freq_2dif={}
for barcode,similar_barcodes in dic_2dif.items():
  barcode_freq = barcodes_dic[barcode]
  list_freq_values = [barcodes_dic[similar_bc] for similar_bc in similar_barcodes] 
  if barcode_freq <= max(list_freq_values): 
    try:
      dic_freq_2dif[barcode_freq] +=1
    except:
      dic_freq_2dif[barcode_freq] =1
dic_freq_2dif_list = sorted(list(dic_freq_2dif.items())) #We order dic_freq_2dif by frequencies:
dic_freq_2dif = dict(dic_freq_2dif_list)
print("\t\tFrequencies of barcodes with 2 differences from other barcodes with equal or higher frequency: ", dic_freq_2dif)

print("\tStep 4 done\n")


print("\tStep 5 => Creating dict with frequencies of the barcodes")
dict_freq_barcodes = defaultdict(int)
for value in barcodes_dic.values():
    dict_freq_barcodes[value] += 1
print("\tStep 5 done\n")


print("\tStep 6 => Plotting frequencies of barcodes with 1 or 2 differences from other barcodes with equal or higher frequency")
plot_file  = out_prefix + "_differences.png"
dic_freq_1dif_percent = {frequence: dic_freq_1dif[frequence] / dict_freq_barcodes[frequence] * 100 for frequence in dic_freq_1dif.keys()}
dic_freq_2dif_percent = {frequence: dic_freq_2dif[frequence] / dict_freq_barcodes[frequence] * 100 for frequence in dic_freq_2dif.keys()}
fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.scatter(list(dict_freq_barcodes.keys()), list(dict_freq_barcodes.values()), marker='o', label='Frequency of barcodes')# Add the data from each dictionary as a scatter plot with circles, squares, and triangles
ax2.plot(list(dic_freq_1dif_percent.keys()), list(dic_freq_1dif_percent.values()), 's-', label='% bc with 1-nt difference from bc >=freq', color='red')
ax2.plot(list(dic_freq_2dif_percent.keys()), list(dic_freq_2dif_percent.values()), '^-', label='% bc with 2-nt difference from bc >=freq', color='green')
ax.set_xlabel('Raw sequences per barcode')
ax.set_ylabel('Number of barcodes')
ax2.set_ylabel('%')
ax.set_title('Assesment of offspring barcodes')
plt.xscale('log')
ax2.set_ylim(0, 100.1)
ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.6))
ax2.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))
plt.savefig(plot_file, bbox_inches='tight')
plt.close()
print("\tStep 6 done\n")


print("\n")
print("######################  END  ######################")