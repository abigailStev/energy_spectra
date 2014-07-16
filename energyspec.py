import argparse
import numpy as np
import itertools

"""
		energyspec.py

Takes CCF amplitudes+mean energy spectrum of a specific time bin and writes them
to a file.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""
###############################################################################
def output(in_file, out_file, bin_num, filt_ccf_amps):
	"""
			output
	
	Writes the energy spectrum to a file. This output file is then used as input 
	for the FTOOLS script ascii2pha.
	
	Passed: in_file - Name of input file with CCF amplitudes per energy bin.
			out_file - Name of output file to write energy spectrum to.
			bin_num - The CCF phase bin to get energy spectrum for.
			filt_ccf_amps - CCF amplitude per energy bin, filtered in frequency.
			
	Returns: nothing
			
	"""
	
	print "Output sent to %s" % out_file

	with open(out_file, 'w') as out:
		for item in filt_ccf_amps:
			out.write("%.5f\t%.8f\n" % (item, np.absolute(item*0.1)))

	## End of function 'output'


###############################################################################
def main(in_file, out_file, bin_num, filt_loc):
	""" 
			main
			
	Finds the time bin desired, gets the filtered CCF amplitudes, sends to 
	output.
	
	Passed: in_file_list - Name of input file with CCF amplitudes per energy 
				bin.
			out_file - Name of output file to write energy spectrum to.
			bin_num - The CCF phase bin to get energy spectrum for.
			filt_loc - 1 if filtered ccf is in first half of columns, 2 if in 
				second half of columns in 'in_file'.

	Returns: nothing
	
	"""
	assert bin_num >= 0
# 	bin_num = -1
	ccf_amps = np.zeros(128, dtype=float)

	with open(in_file, 'r') as f:
		for line in f:
			if line[0].strip() != "#":
# 				print line
				line = line.strip().split()
				if int(line[0]) == bin_num:
# 					print "Bin found!"
					for i in range (0, 128):
						ccf_amps[i] = float(line[i+1])
					break
		else:
			print "\n\tERROR: Phase bin not found. Check that it is within the range of the file. Exiting."
			exit()
	## End of with-block
	
	
	if filt_loc == 1:
		filt_ccf_amps = ccf_amps[0:64]
	else:
		filt_ccf_amps = ccf_amps[64:128]

	
	
# 	print np.shape(filt_ccf_amps)
	
	output(in_file, out_file, bin_num, filt_ccf_amps)

	## End of function 'main'


###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Takes CCF amplitudes + mean \
		energy spectrum of a specific time bin and writes them to a file.')
	parser.add_argument('-i', '--infile', required=True, dest='infile', 
		help="The full path of the (ASCII/txt/dat) input file listing the CCF \
		amplitudes per energy bin.")
	parser.add_argument('-o', '--outfile', required=True,  dest='outfile', 
		help="The full path of the (ASCII/txt/dat) output file to write the \
		cross-correlation function to.")
	parser.add_argument('-b', '--bin', required=True, type=int, dest='bin_num',  
		help="The phase bin number of the CCF to select energy bins for.")
	parser.add_argument('-f', '--filt_loc', type=int, default=1, \
		choices=range(1, 3), dest='filt_loc', help="1 if filtered CCF is in \
		first half of columns, 2 if in second half of columns in the CCF \
		amplitude table.")
	args = parser.parse_args()

	main(args.infile, args.outfile, args.bin_num, args.filt_loc)

## End of program 'energyspec.py'
