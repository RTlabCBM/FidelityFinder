# encoding: UTF-8


# Imports
import sys, os, math
import json
import seaborn as sbn
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from Bio import SeqIO
from collections import defaultdict
from scipy import stats
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.AlignInfo import SummaryInfo

########################################################################################
print()
print("\t", "*" * 75)
print()
print("\tPython script to calculate consensus sequence by barcode")
print()
print("\t-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .")
print("\t  \   \ /   \   \ /   \   \ /   \   \ /   \   \ /   \   \ /   \   \ / ")
print("\t / \   \   / \   \   / \   \   / \   \   / \   \   / \   \   / \   \  ")
print("\t~   `-~ `-`   `-~ `-`   `-~ `-`   `-~ `-`   `-~ `-`   `-~ `-`   `-~ `-")
print()
print("\t", "*" * 75, "\n")
########################################################################################

__doc__ = """
SYNOPSIS

 Python script to calculate consensus sequences by barcode 

    consensus_barcode.py -f <file> -r <file> -p <str> -q <str> -b <int> -i <int> -s <int> -t <float> [-o,--output] [-h,--help] [-v,--version]

DESCRIPTION

 Parameters:
    -f,--fasta		            fasta file [required]
    -r,--reference              reference sequence (each position of the barcode must be indicated with "N") [required]
    -p,--fw_primer	            forward primer [default='CAGGAGCCGATAGACAAGGAAC']
    -q,--rv_primer	            reverse primer [default='GGAATGGATGGCCCAAAAGTTAAACTG']
    -s,--superior	            upper cut-off of sequences per barcode [default=100000]
    -i,--inferior	            lower cut-off of sequences per barcode [default=4]
    -t,--threshold_consensus    threshold to construct consensus sequences [default=0.9]
    -o,--output_prefix		    output prefix [default='consensus']
    -b,--bc_size	            barcode size (bp) [default=14]
    -v,--version	            Shows script version
    -h,--help		            Shows help options


AUTHORS

    This is a script modified by Javier Martinez del Río (javier.martinez@cbm.csic.es) and 
    Estrella Frutos-Beltrán (efrutos@cbm.csic.es) from the version 1 of the script created 
    by the following authors:
    
    Genomics & NGS Facility (CBMSO-CSIC) | Eva Castillo    <ecastillo@cbm.csic.es>
                                         | Eva Sacristán  <esacristan@cbm.csic.es>
                                         | Sandra González  <sandra.g@cbm.csic.es>
                                         | Ramón Peiró-Pastor <rpeiro@cbm.csic.es>
                    
LICENSE

    Copyright (c) 2018 CBMSO
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
__author__ = "Eva Castillo, Eva Sacristán, Sandra González, Ramón Peiró-Pastor, Javier Martínez del Río, Estrella Frutos-Beltrán"
__version__ = 'v2.0.0'

## Functions
##----------

def help():
    print(globals()['__doc__'])
    sys.exit(1)

def paramerror(param, message):
    """Gives an error message if argument is not valid"""
    print("\n\t########## ERROR ########## \n'{} {}'\n".format(param, message))
    sys.exit(1)

def send_error_message(message):
    """Gives an error message if argument is not valid"""    
    print("\n\t########## ERROR ########## \n{}\n".format(message))
    sys.exit(1)

def version():
    """Print the version of this script"""
    print("\n{}".format(script_name))
    print("Version: {}".format(__version__))
    print("Authors: {}\n".format(__author__))
    sys.exit(1)

def find_barcode(reference_sequence, query_sequence):
    # reference_sequence is the insert of the library sequenced indicating each letter of the barcode zone as "N"
    # query_sequence is the raw sequence of the insert obtained by NGS 
    # Align reference and query sequences:
    mafft_cline = MafftCommandline("mafft", input="-", auto=True, thread=30)
    stdout, stderr = mafft_cline(stdin=f">reference_sequence\n{reference_sequence}\n>query_sequence\n{query_sequence}")
    align = AlignIO.read(StringIO(stdout), "fasta")
    ref_aligned = str(align[0].seq)
    query_aligned = str(align[1].seq)
    # Select barcode sequence (it corresponds to the part of query_sequence that aligns with the letters "N" of reference_sequence)
    barcode = ""
    for ref_nuc, query_nuc in zip(ref_aligned, query_aligned):
        if ref_nuc == "n" and  query_nuc != "-":
            barcode += query_nuc
    return barcode


## Parameters catching
##--------------------

# defaults
bc_size    = 14
inferior   = 4
superior   = 1000000
out_prefix = 'undefined'
fasta_file = ''
fprimer    = 'CAGGAGCCGATAGACAAGGAAC'
rprimer    = 'GGAATGGATGGCCCAAAAGTTAAACTG'
threshold_consensus = 0.9


print("\tStep 1 => Parameters catching")
script_name = sys.argv.pop(0)
if len(sys.argv) == 0:
  send_error_message("You have to provide some parameters !!!")
  version()
else:
  while len(sys.argv) > 0:
        param = sys.argv.pop(0)
        if param == '-f' or param == '--fasta':
            fasta_file = sys.argv.pop(0)
            if not os.path.exists(fasta_file):
                paramerror(fasta_file, "file does not exist !!!")
        elif param == '-r' or param == '--reference':
            reference_file = sys.argv.pop(0)
            if not os.path.exists(reference_file):
                paramerror(reference_file, "file does not exist !!!")                
        elif param == '-b' or param == '--bc_size':
            bc_size = int(sys.argv.pop(0))
        elif param == '-i' or param == '--inferior':
            inferior = int(sys.argv.pop(0))
        elif param == '-s' or param == '--superior':
            superior = int(sys.argv.pop(0))
        elif param == '-o' or param == '--output_prefix':
            out_prefix = sys.argv.pop(0)
        elif param == '-p' or  param == '--fw_primer':
            fprimer = sys.argv.pop(0)
            fprimer = fprimer.upper()	# Preventing upper and lower case mixture
        elif param == '-q' or  param == '--rv_primer':
            rprimer = sys.argv.pop(0)
            rprimer = rprimer.upper()	# Preventing upper and lower case mixture
        elif param == '-t' or  param == '--threshold_consensus':
            threshold_consensus = float(sys.argv.pop(0))
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

## Main program
##-------------
print("\tStep 2 => Open & read FASTA file | Barcode search")
fasta_hand = SeqIO.parse(fasta_file, "fasta")
reference_sequence = next(SeqIO.parse(reference_file, "fasta")).seq
fasta_dic = defaultdict(list)
bc_lengths_dic = defaultdict(int)
accepted_sequences = 0
rejected_sequences_bc_size = 0
rejected_sequences_no_align = 0
for record in fasta_hand:
  read     = str(record.seq)
  try:
    barcode  = find_barcode(reference_sequence, read)
    bc_lengths_dic[len(barcode)] += 1
    if len(barcode) == bc_size:
        fasta_dic[barcode].append(read)
        accepted_sequences += 1
    else:
        rejected_sequences_bc_size += 1
  except Exception as e:
    rejected_sequences_no_align += 1
    print(e)

print("\t\tTotal initial sequences: ", accepted_sequences + rejected_sequences_bc_size + rejected_sequences_no_align)
print("\t\tAccepted sequences: ", accepted_sequences)
print("\t\tRejected sequences (different barcode length): ", rejected_sequences_bc_size)
print("\t\tRejected sequences (no alignment to the reference sequence): ", rejected_sequences_no_align)
print("\t\tBarcodes lengths and their frequency: ", dict(bc_lengths_dic))
print("\t\tTotal unique barcodes with expected length found: ", len(fasta_dic.keys()))

if accepted_sequences < 1:
    send_error_message(f'No sequences could be found (sample: {out_prefix}), please check the input parameters of {script_name}\nHave you chosen the right input parameters?')

print("\tStep 2 done\n")





print("\tStep 3 => Calculating consensus sequences")

# files
fasta_tmp = "tmp.fa"
xls_file                = out_prefix + "_consensus.xls"
discarded_file          = out_prefix + "_discarded.txt"
profiles_file           = out_prefix + "_consensus.prf"
consensus_file          = out_prefix + "_consensus.fna"
barcodes_file           = out_prefix + "_barcodes.json"

# handlers
xls_hand                = open(xls_file, "w")
discarded_hand          = open(discarded_file, "w")
prf_hand                = open(profiles_file, "w")
consensus_hand          = open(consensus_file, "w")
barcodes_hand           = open(barcodes_file, 'w')

# header for xls file
header_string = "#barcode\tconsensus\tnum_seqs\n"
xls_hand.write(header_string)

# store info to draw plot
plot_dic  = {}
plot2_dic = {}
barcodes_dic = {}

for barcode in fasta_dic.keys():
  counter  = 0
  seqList  = fasta_dic[barcode]
  nseqs    = len(seqList)
  barcodes_dic[barcode] = nseqs
  try:
    plot_dic[nseqs] += 1
  except KeyError:
    plot_dic[nseqs] = 1
  if nseqs>inferior and nseqs<superior:
    try:
      plot2_dic[nseqs] += 1
    except KeyError:
      plot2_dic[nseqs] = 1
    tmp_hand = open(fasta_tmp, "w")
    for seq in seqList:
      string = ">" + barcode + "_" + str(counter) + "\n" + seq + "\n"
      tmp_hand.write(string)
      counter += 1
    tmp_hand.close()
    mafft_cline = MafftCommandline("mafft", input=fasta_tmp, thread=30)
    stdout, stderror= mafft_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    align_sum = SummaryInfo(align)
    consensus = str(SummaryInfo.gap_consensus(align_sum, threshold=threshold_consensus, ambiguous='N')).upper()
    matrix = str(SummaryInfo.pos_specific_score_matrix(align_sum)).upper()
    os.remove(fasta_tmp)
    string = ">" + barcode + " nseqs=" + str(nseqs) + "\n" + consensus + "\n"
    consensus_hand.write(string)
    string = "barcode=" + barcode + " nseqs=" + str(nseqs) + "\n" + matrix + "\n"
    prf_hand.write(string)
    string = barcode + "\t" + consensus + "\t" + str(nseqs) + "\n"
    xls_hand.write(string)
  else:
    string = "barcode=" + barcode + " nseqs=" + str(nseqs) + "\n"
    discarded_hand.write(string)
    
json.dump(barcodes_dic, barcodes_hand)
barcodes_hand.close()
consensus_hand.close()
prf_hand.close()
xls_hand.close()
discarded_hand.close()

print("\t\tSequences per barcode and their respective frequencies: ", plot_dic)
print(f"\t\tTotal number of consensus sequences constructed (using {inferior} as inferior cutoff value): ", sum(plot2_dic.values())) 
if len(plot2_dic.keys()) < 1:
    send_error_message(f'No consensus sequences could be constructed using {inferior} as inferior cutoff value')
    
print("\tStep 3 done\n")



print("\tStep 4 => Drawing plot")
plot_file  =  out_prefix + "_consensus.png"
x = plot_dic.keys()
y = plot_dic.values()
plt.scatter(x, y, marker="o")
plt.yscale('log')
titlename = 'Primer ID distribution for ' + out_prefix + ' reads (no cut-off)'
plt.title(titlename)
plt.xlabel('Raw Sequence Reads per Unique Primer ID')
plt.ylabel('Number of Distinct Primer IDs')
plt.savefig(plot_file)
plt.close()

plot2_file  = out_prefix + "_cutoff_consensus.png"
x = plot2_dic.keys()
y = plot2_dic.values()
plt.scatter(x, y, marker="o")
plt.yscale('log')
titlename = 'Primer ID distribution for ' + out_prefix + ' reads cut-off (' + str(inferior) + ", " + str(superior) + ")"
plt.title(titlename)
plt.xlabel('Raw Sequence Reads per Unique Primer ID')
plt.ylabel('Number of Distinct Primer IDs')
plt.savefig(plot2_file)
plt.close()

print("\tStep 4 done\n")

print("\tJOB DONE !!! Have a nice & pythonic day!!! (^_^)\n")