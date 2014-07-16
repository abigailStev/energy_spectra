import argparse
import numpy as np
import itertools
from astropy.io import fits

import tools

"""
		mean_spectrum

Computes the mean spectrum from many event.pha files and adds it to the ccf 
amplitude.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""
###############################################################################
def output(out_file, bins, out_tab):
	"""
			output
	
	
	
	Passed: out_file
			bins
			out_tab
			
	Returns: nothing
			
	"""
	
	print "Output sent to %s" % out_file

	with open(out_file, 'w') as out:
		out.write("# ")
		for j in xrange(0, len(bins)):
			out.write("\n%d" % bins[j])
			for i in xrange(0, 64):
				out.write("\t%.6f" % out_tab[j][i])
	# 		for i in xrange(0, 64):
	# 			out.write("\t%.5f" % ccf_error[j][i].real)
	
		## End of for-loops
	## End of with-block
			

## End of function 'output'


###############################################################################
def main(in_file_list, ccf_file, out_file):
	""" 
			main
			
	
			
	Passed: in_file_list
			ccf_file
			out_file

	Returns: nothing
	
	"""

	data_files = [line.strip() for line in open(in_file_list)]
	num_files = len(data_files)
# 	print num_files
	total_counts = np.zeros(64, dtype=np.float64)
	total_counts_err = np.zeros(64, dtype=np.float64)
	total_exposure = 0
	
	for in_file in data_files:
		file_hdu = fits.open(in_file)
		fits_data = file_hdu[1].data
		counts = fits_data.field('COUNTS')
		counts_err = fits_data.field('STAT_ERR')
		exposure = float(file_hdu[1].header['EXPOSURE'])
		total_exposure += exposure
# 		print exposure
		file_hdu.close()
		
		total_counts += counts
		total_counts_err += counts_err
	## End of for-loop
	
# 	total_countrate /= float(num_files)
# 	total_countrate_err /= float(num_files)

	total_countrate = total_counts / total_exposure
	total_countrate_err = total_counts_err / total_exposure
	print "Total countrate:", total_countrate
# 	print total_countrate_err

	ccf_table = np.loadtxt(ccf_file, comments='#')
	print "Shape of ccf table:", np.shape(ccf_table)
	
	bins = ccf_table[:,0].astype(int)
	ccf_amps = ccf_table[:,1:65]
	ccf_errs = ccf_table[:,66:]
	print "Shape of ccf amps:", np.shape(ccf_amps)
# 	print "CCF elt:", ccf_amps[1,0]
	
	out_tab = np.add(ccf_amps, total_countrate)
	print "Shape of mean spectra + ccf:", np.shape(out_tab)
# 	print "newtemp elt:", new_temp[1,0]
	
	
	output(out_file, bins, out_tab)

## End of function 'main'


###############################################################################
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Computes the mean spectrum from many event.pha files and adds it to the ccf amplitude.')
	parser.add_argument('-i', '--infile_list', required=True, 
		dest='infile_list', help='File containing list of event.pha input \
		files. One file per line.')
	parser.add_argument('-o', '--outfile', required=True, dest='outfile', 
		help='')
	parser.add_argument('-c', '--CCFfile', required=True, dest='ccffile', 
		help='')
	args = parser.parse_args()

	main(args.infile_list, args.ccffile, args.outfile)

## End of program 'mean_spectrum.py'
