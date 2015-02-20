import argparse
import numpy as np
import itertools
from astropy.io import fits
from tools import read_obs_time  # https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014-2015"
__description__ = "Takes CCF amplitudes+mean energy spectrum of a specific time\
 bin and writes them to a file. Can indicate whether to produce spectra that \
are mean+ccf, ccf, or mean."


"""
		energyspec.py

Written in Python 2.7.

"""

################################################################################
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
	
# 	print "Output sent to %s" % out_file
	
	with open(out_file, 'w') as out:
		for i in xrange(0, 64):
			out.write("%d\t%.6e\t%.6e\n" % (i, amps[i], err[i]))
# 			out.write("%d\t%.6e\t%.6e\n" % (i, amps[i], amps[i]*.1))
	## End of with-block
	
## End of function 'output'


################################################################################
def get_mean_count_rate(string):
	"""
			get_mean_count_rate
	
	Passed: string - The comment line from the ccf file with mean rate per 
				energy channel.
	
	Returns: an array of the mean rates per energy channel
	
	"""
	
# 	print string
	start_index = string.index('[')
# 	print start_index
# 	print string[start_index+1:-1]
	return np.asarray(string[start_index+1:-1].split(', '), dtype=np.float64)

## End of 'get_mean_count_rate'	


################################################################################
def main(in_file, out_file, bin_num, spec_type):
	""" 
			main
			
	Finds the time bin desired, gets the filtered CCF amplitudes, sends to 
	output.
	
	"""
	assert bin_num >= 0, "ERROR: Bin number must be >= 0."
# 	bin_num = -1
	ccf_amps_and_err = np.zeros(128, dtype=np.float64)
	mean_count_rate = np.zeros(64, dtype=np.int32)
	
	if in_file[-3:].lower() == 'dat':
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
				raise Exception("ERROR: Phase bin not found. Check that it is within the range of the file.")
		## End of with-block
	
		ccf_amps = ccf_amps_and_err[0:64]
		ccf_err = ccf_amps_and_err[64:128]
		obs_time = read_obs_time(in_file)
	
	elif in_file[-4:].lower() == 'fits':
		file_hdu = fits.open(in_file)
		table = file_hdu[1].data
		obs_time = file_hdu[0].header['EXPOSURE']
		mean_count_rate = get_mean_count_rate(file_hdu[0].header['RATE_CI'])
		mean_count_rate[np.where(mean_count_rate < 0.0)] = 0
		file_hdu.close()
		
		time_bin_mask = table.field('TIME_BIN') == bin_num
		table_i = table[time_bin_mask]
		ccf_amps = table_i.field('CCF')
		ccf_err = table_i.field('ERROR')
	else:
		raise Exception("ERROR: Input file must have extension .dat or .fits.")
	## End of if/elif file type
	
	mean_err = np.sqrt(mean_count_rate * obs_time) / obs_time
	amps = []
	err = []
# 	print spec_type
	if spec_type == 0:
		amps = np.add(ccf_amps, mean_count_rate) 
		err = np.sqrt(np.add(np.square(ccf_err), np.square(mean_err)))
	elif spec_type == 1:
		amps = ccf_amps
		err = ccf_err
	elif spec_type == 2:
		amps = mean_count_rate
		err = mean_err
	else:
		raise Exception("ERROR: Spectrum type not a valid option.")
	## End of if/else spectrum type

	output(out_file, bin_num, amps, err)

## End of function 'main'


################################################################################
if __name__ == "__main__":
	
	##############################################
	## Parsing input arguments and calling 'main'
	##############################################
	
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
	parser.add_argument('-s', '--spec', type=int, dest='spec_type', 
		choices=range(0, 3), help="Indicating the type of spectrum to produce. \
		0 for mean+ccf, 1 for ccf, 2 for mean.")
	args = parser.parse_args()

	main(args.infile, args.outfile, args.bin_num, args.spec_type)

## End of program 'energyspec.py'
################################################################################
