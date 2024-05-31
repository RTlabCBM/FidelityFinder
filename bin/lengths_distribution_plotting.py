# encoding: UTF-8


__doc__ = """
SYNOPSIS

  Python script to create bar plots depicting the distribution of sequence lengths. 

    lengths_distribution_plotting.py -i <file> -o <str>  [-v,--version] [-h,--help]        

DESCRIPTION

 This script reads sequence length data from an input file and generates two bar plots:
 1. A standard bar plot with sequence lengths on the x-axis and corresponding frequencies on the y-axis.
 2. A log-scaled bar plot for better visualization of frequency distribution.

 Parameters:
    -i,--input_file         file [required]
    -o,--output_prefix		output prefix [default='undefined']
    -v,--version	        Shows script version
    -h,--help		        Shows help options
    
 Input file:
 Input file must contain rows where the first element of each row is the length of a sequence and the second element is the frequency of that length. Example:
    450 2
    451;16
    452,7
    
  Output files:  
  Two plot files are generated with the specified output prefix:
    1. {output_prefix}_lengths_distribution.png
    2. {output_prefix}_lengths_distribution_logscale.png


AUTHORS

    Javier Martinez del Río (javier.martinez@cbm.csic.es; javier.mardelrio@gmail.com)   
    
"""

__version__ = 'v1.0.0'


## Imports
##----------
import argparse
import sys, os
import re
import matplotlib.pyplot as plt
plt.switch_backend('agg')



## Functions
##----------

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input_file', required=True, help='Input file')
    parser.add_argument('-o', '--output_prefix', default="undefined", help='Output prefix')
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    return parser.parse_args()
    
def send_error_message(message):
    """Gives an error message"""    
    print("\n\t########## ERROR ########## \n{}\n".format(message))
    sys.exit(1)
    
def read_data(input_file):
    with open(input_file, 'r') as file:
        content = file.read()
    rows = content.strip().split('\n')
    lengths = []
    frequencies = []
    for row in rows:
        numbers = re.split(r'\D+', row)
        lengths.append(int(numbers[0]))
        frequencies.append(int(numbers[1]))
    return lengths, frequencies
    
def create_plot(lengths, frequencies, output_file, log_scale=False):
    fig = plt.figure(figsize=(10, 6))
    plt.bar(lengths, frequencies, width=2)
    plt.xlabel('Sequence length (nt)')
    plt.ylabel('Number of sequences')
    title = 'Distribution of sequence lengths'
    if log_scale:
        title += ' (log scale)'
        plt.yscale('log')
    plt.title(title)
    plt.savefig(output_file)
    plt.close()


    
## Main program
##-------------

def main():

    ## Step 1: Parameters catching
    ##-------------------- 

    print("\tStep 1 => Parameters catching")

    args = parse_arguments()
    input_file = args.input_file
    output_prefix = args.output_prefix
    
    if not os.path.isfile(input_file):
        send_error_message(f'Input file "{input_file}" of {__file__} does not exist. Please provide a valid file path.')
        
    print("\tStep 1 done\n")


    ## Step 2: Reading input file
    ##-------------------- 
    
    print("\tStep 2 Read input file\n")
    
    lengths, frequencies = read_data(input_file)

    print("\tStep 2 done\n")
    
    
    ## Step 3: Creating graphs
    ##--------------------  
    
    print("\tStep 3 Creating graphs\n")
 
    # FIRST GRAPH (standard axis)
    plot_file = f"{output_prefix}_lengths_distribution.png"
    create_plot(lengths, frequencies, plot_file, log_scale=False)

    # SECOND GRAPH (log axis)
    plot_file_log = f"{output_prefix}_lengths_distribution_logscale.png"
    create_plot(lengths, frequencies, plot_file_log, log_scale=True)   
    
    print("\tStep 3 done\n")
    

if __name__ == "__main__":
    main()