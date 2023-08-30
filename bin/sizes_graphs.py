# encoding: UTF-8


# Imports
import matplotlib.pyplot as plt
import numpy as np
import sys, os, math


__doc__ = """
SYNOPSIS

 Python script to create plots with lengths of sequences. 
 

    sizes_graphs.py -i <file> -o <str> [-h,--help]        

DESCRIPTION

 Parameters:
    -i,--input_file         file [required]
    -o,--output_prefix		output prefix [default='']
    -h,--help		        Shows help options

 Input file must contain rows where the first element of each row is the length of a sequence and the second element is the frequency of that length. 
 For example:
    450 2
    451 16
    452 7
    
AUTHORS

    Javier Martinez del Río (javier.martinez@cbm.csic.es) and Estrella Frutos Beltrán (efrutos@cbm.csic.es) 
    
"""


## Functions
##----------


def help():
    print(globals()['__doc__'])
    sys.exit(1)

def paramerror(param, message):
  """Gives an error message if argument is not valid"""
  print("\tERROR: '", param, "'", message, "\n")
  sys.exit(1)
  
def paramerrorB(message):
    """Gives an error message if argument is not valid"""
    print("\tERROR:", message, "\n")
    sys.exit(1)


## Parameters catching
##--------------------


# defaults
output_name = ""

print("\tStep 1 => Parameters catching")
script_name = sys.argv.pop(0)
if len(sys.argv) == 0:
  paramerrorB("You have to provide some parameters !!!")
else:
  while len(sys.argv) > 0:
        param = sys.argv.pop(0)
        if param == '-i' or param == '--input_file':
            input_file = sys.argv.pop(0)
            if not os.path.exists(input_file):
                paramerror(input_file, "input file does not exist !!!")
        elif param == '-o' or param == '--output_prefix':
            output_prefix = sys.argv.pop(0)
        elif param == '-h' or param == '--help':
            help()
            sys.exit(1)
        else:
            if param == script_name:
                continue
            else:
                paramerror(param, "parameter not recognized !!!")
print("\tStep 1 done\n")


## Main program
##--------------------

print("\tStep 2 Creating graphs\n")

#FIRST GRAPH (standard axis)
#For each row of the input file, we select the first value (length of the sequence) and the second value (frequency of the length)
#Then we graph length of sequences (x axis) vs frequency of the lengts (y axis)
plot_file  = output_prefix + "_sequences_sizes_after_merging.png"
with open(input_file, 'r') as file:
    content = file.read()
rows = content.strip().split('\n')
x = []
y = []
for row in rows:
    x.append(int(row.split()[0]))
    y.append(int(row.split()[1]))    
fig = plt.figure(figsize=(10, 6))
plt.bar(x,y, width = 2)
plt.xlabel('Size (bp)')
plt.ylabel('Number of sequences')
plt.title('Sequences sizes after merging')
plt.savefig(plot_file)
plt.close()


#SECOND GRAPH (log axis)
#For each row of the input file, we select the first value (length of the sequence) and the second value (frequency of the length)
#Then we graph length of sequences (x axis) vs frequency of the lengts (y axis)
plot_file_log  = output_prefix + "_sequences_sizes_after_merging_logscale.png"   
fig = plt.figure(figsize=(10, 6))
plt.bar(x,y, width = 2)
plt.xlabel('Size (bp)')
plt.ylabel('Number of sequences')
plt.title('Sequences sizes after merging (log scale)')
plt.yscale('log')
plt.savefig(plot_file_log)
plt.close()

print("\tStep 2 done\n")
print("\tEND")
