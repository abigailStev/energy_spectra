import argparse
import numpy as np
import itertools
from tools import read_obs_time

"""
		energyspec.py

Takes CCF amplitudes+mean energy spectrum of a specific time bin and writes them
to a file.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""
###############################################################################
def output(out_file, bin_num, amps, err):
	"""
			output
	
	Writes the energy spectrum to a file. This output file is then used as input 
	for the FTOOLS script ascii2pha.
	
	Passed: out_file - Name of output file to write energy spectrum to.
			bin_num - The CCF phase bin to get energy spectrum for.
			amps - Amplitude per energy bin
			err - Error per energy bin
			
	Returns: nothing
			
	"""
	pass
	
	print "Output sent to %s" % out_file
	
	with open(out_file, 'w') as out:
		for i in xrange(0, 64):
			out.write("%d\t%.6e\t%.6e\n" % (i, amps[i], err[i]))
# 			out.write("%d\t%.6e\t%.6e\n" % (i, amps[i], amps[i]*.1))
	## End of with-block
	
## End of function 'output'


###############################################################################
def get_mean_count_rate(string):
	"""
			get_mean_count_rate
	
	Passed: string - The comment line from the ccf file with mean rate per 
				energy channel.
	
	Returns: an array of the mean rates per energy channel
	
	"""
	pass
	
# 	print string
	start_index = string.index('[')
# 	print start_index
# 	print string[start_index+1:-1]
	return np.asarray(string[start_index+1:-1].split(', '), dtype=np.float64)

## End of 'get_mean_count_rate'	


###############################################################################
def main(in_file, out_file, bin_num):
	""" 
			main
			
	Finds the time bin desired, gets the filtered CCF amplitudes, sends to 
	output.
	
	Passed: in_file - Name of input file with CCF amplitudes per energy bin.
			out_file - Name of output file to write energy spectrum to.
			bin_num - The CCF phase bin to get energy spectrum for.

	Returns: nothing
	
	"""
	pass
	assert bin_num >= 0
# 	bin_num = -1
	ccf_amps_and_err = np.zeros(128, dtype=np.float64)
	mean_count_rate = np.zeros(64, dtype=np.int32)

	with open(in_file, 'r') as f:
		for line in f:	
			if line[0].strip() != "#":
# 				print line
				line = line.strip().split()
				if int(line[0]) == bin_num:
# 					print "Bin found!"
					for i in xrange(0, 128):
						ccf_amps_and_err[i] = float(line[i+1])
					break
			else:
				if "Mean" in line.strip() and \
					"count rate" in line.strip() and \
					"ci" in line.strip():
					mean_count_rate = get_mean_count_rate(line.strip())
		else:
			print "\n\tERROR: Phase bin not found. Check that it is within the range of the file. Exiting."
			exit()
	## End of with-block
	
	ccf_amps = ccf_amps_and_err[0:64]
# 	print ccf_amps
	ccf_err = ccf_amps_and_err[64:128]
	
	obs_time = read_obs_time(in_file)
	mean_err = np.sqrt(mean_count_rate * obs_time) / obs_time
# 	print np.shape(ccf_amps)
# 	print ccf_amps

# 	amps = np.add(ccf_amps, mean_count_rate) 
	amps = mean_count_rate
# 	print mean_count_rate[20:23]
	
# 	err = np.sqrt(np.add(np.square(ccf_err), np.square(mean_err)))
	err = mean_err

	output(out_file, bin_num, amps, err)

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
	args = parser.parse_args()

	main(args.infile, args.outfile, args.bin_num)

## End of program 'energyspec.py'
