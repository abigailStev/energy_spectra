import argparse
import numpy as np
from astropy.io import fits

#########################################################################################
def make_mean_spectrum(obsID_list, dt):
	"""
	Using rules for combining errors to get the countrate error.
	
	if T = |k| * (A + B + ...)
	then 
	E_T = |k| * sqrt( E_A**2 + E_B**2 + ... )
	"""
	mean_countrate = np.zeros(64)
	mean_error = np.zeros(64)
	obsIDs = [line.strip() for line in open(obsID_list)]
	num_files = float(len(obsIDs))
	
	for obsID in obsIDs:
		file = "/Users/abigailstevens/Reduced_data/P70080/"+str(obsID)+"/event.pha"
		fits_hdu = fits.open(file)
		header = fits_hdu[1].header	
		data = fits_hdu[1].data
		fits_hdu.close()
		# print header.keys()
		chan = data.field(0)
		countrate = data.field(1) * dt
		error = data.field(2) * dt
		print countrate
		mean_countrate += countrate
		mean_error += error ** 2
		
	mean_countrate /= num_files
	mean_error = np.sqrt(error) / num_files
	
	return mean_countrate, mean_error


#########################################################################################
def main(phase_file, obsID_list, out_file):
	"""
	if T = A / B
	then
	E_T = T * sqrt( (E_A / A)**2 + (E_B / B)**2 + ...)
	
	
	phase_file = "Dropbox/Research/power_spectra/out_es/140610_t1_4sec_pbin_25.dat"
 	mean_file = "Reduced_data/P70080/70080-01-02-01/event.pha"
	out_file = "Dropbox/Research/power_spectra/out_es/140610_t1_4sec_pbin_25_ratio.dat"
	"""
	pass

	dt = 1./8192.0
	old_settings = np.seterr(all='ignore')

	print "Phase spectrum: %s" % phase_file

	phase_table = np.loadtxt(phase_file, comments='#')
	phase_chan = np.arange(64, dtype='int8')
	phase_countrate = phase_table[:,0]
	phase_error = phase_table[:,1]

	mean_countrate, mean_error = make_mean_spectrum(obsID_list, dt)

	countrate_ratio = np.nan_to_num(phase_countrate / mean_countrate)
	error_ratio = (phase_error / phase_countrate) + (mean_error / mean_countrate)
	error_ratio = np.nan_to_num(countrate_ratio * np.sqrt(error_ratio))

	ratio_spectrum = np.column_stack((countrate_ratio, error_ratio))
# 	print np.shape(ratio_spectrum)

	print "Out file: %s" % out_file

	np.savetxt(out_file, ratio_spectrum, fmt='%.8e\t%.8e')
	

#########################################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('phase_file', help="")
	parser.add_argument('obsID_list', help="")
	parser.add_argument('out_file', help="")
	args = parser.parse_args()

	main(args.phase_file, args.obsID_list, args.out_file)