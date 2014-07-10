import argparse
import math
import numpy as np
from scipy import fftpack
import itertools

"""
		energyspec.py

Takes CCF amplitudes of a specific time bin and writes them to a file.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

"""
#########################################################################################
def output(in_file, out_file, bin_num, filt_ccf_amps):
	"""
			output
	
	Writes the energy spectrum to a file. This output file is then used as input for 
	the FTOOLS script ascii2pha.
	
	Passed: in_file - Name of input file with CCF amplitudes per energy bin.
			out_file - Name of output file to write energy spectrum to.
			bin_num - The CCF phase bin to get energy spectrum for.
			filt_ccf_amps - CCF amplitude per energy bin, filtered in frequency.
			
	Returns: nothing
			
	"""
	
	print "Output sent to %s" % out_file

	with open(out_file, 'w') as out:
		for item in filt_ccf_amps:
			out.write("%.5f\t%.8f\n" % (item, math.fabs(item*0.1)))

	## End of function 'output'


#########################################################################################
def main(in_file, out_file, bin_num, filt_loc):
	""" 
			main
			
	Finds the time bin desired, gets the filtered CCF amplitudes, sends to output.
	
	Passed: in_file_list - Name of input file with CCF amplitudes per energy bin.
			out_file - Name of output file to write energy spectrum to.
			bin_num - The CCF phase bin to get energy spectrum for.
			filt_loc - 1 if filtered ccf is in first half of columns, 2 if in second half
				of columns in 'in_file'.

	Returns: nothing
	
	"""
	assert bin_num >= 0
	assert filt_loc == 1 or filt_loc == 2

	ccf_amps = np.zeros(128, dtype=float)

	## Reading only the first line of data to get the start time of the file
	with open(in_file, 'r') as f:
		for line in f:
			if line[0].strip() != "#":
# 				print line
				line = line.strip().split()
				if int(line[0]) == bin_num:
# 					print "Bin found!"
					for i in range (0,128):
						ccf_amps[i] = float(line[i+1])
					break
	
	if filt_loc == 1:
		filt_ccf_amps = ccf_amps[0:64]
	else:
		filt_ccf_amps = ccf_amps[64:128]
	
# 	print np.shape(filt_ccf_amps)
	
	output(in_file, out_file, bin_num, filt_ccf_amps)

	## End of function 'main'


#########################################################################################
if __name__ == "__main__":
	"""
	
	Parsing cmd-line arguments and calling 'main'
	
	"""
	
	parser = argparse.ArgumentParser()
	parser.add_argument('in_file', help="The full path of the (ASCII/txt/dat) input \
		file listing the CCF amplitudes per energy bin.")
	parser.add_argument('out_file', help="The full path of the (ASCII/txt/dat) output \
		file to write the cross-correlation function to.")
	parser.add_argument('bin_num', type=int, help="The phase bin number of the CCF to \
		gather energy bins for.")
	parser.add_argument('filt_loc', type=int, help="1 if filtered ccf is in first half \
		of columns, 2 if in second half of columns in the ccf amplitude table.")
	args = parser.parse_args()

	main(args.in_file, args.out_file, args.bin_num, args.filt_loc)

## End of program 'energyspec.py'
