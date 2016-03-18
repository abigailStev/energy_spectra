#!/bin/bash

################################################################################
## 
## Bash script to run energyspec.py and XSPEC.
##
## Example call: ./run_energyspec.sh J1808 4 16 0 150131
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries) 
## 		  must be installed in order to run this script. 
##
## Written by Abigail Stevens <A.L.Stevens at uva.nl>, 2014-2016
##
################################################################################

## Checking the number of input arguments
if (( $# != 5 )); then
    echo -e "\tUsage: ./run_energyspec.sh <prefix> <dt multiple> <num seconds>"\
            " <testing> <date>\n"
    exit
fi

prefix=$1
dt=$2
numsec=$3
testing=$4
day=$5

################################################################################

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. $HEADAS/headas-init.sh
fi

home_dir=$(ls -d ~)
in_dir="$home_dir/Dropbox/Research/cross_correlation"
exe_dir="$home_dir/Dropbox/Research/energy_spectra"
out_dir="$exe_dir/out_es/${prefix}"
dump_file=dump.txt # Name of dumping file for intermediary steps

data_dir="$home_dir/Reduced_data/$prefix"
# data_dir="$home_dir/Dropbox/Research/sample_data"

bkgd_spec="$data_dir/evt_bkgd_rebinned.pha"
rsp_matrix="$data_dir/PCU2.rsp"
out_root="${prefix}_${day}_t${dt}_${numsec}sec"

ccf_file="$in_dir/out_ccf/${prefix}/${out_root}_adj.fits"
if (( $testing==1 )); then
	ccf_file="$in_dir/out_ccf/${prefix}/test_${out_root}.fits"
fi

tab_ext="dat"
plot_ext="png"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi

########################################
## Reading exposure time of observation
########################################

if [ -e "$ccf_file" ]; then
	if [ "${ccf_file##*.}" == "dat" ]; then
		obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$ccf_file')")
	elif [ "${ccf_file##*.}" == "fits" ]; then
		obs_time=$(python -c "from tools import get_key_val;  print get_key_val('$ccf_file', 1, 'EXPOSURE')")
	fi
	echo -e "EXPOSURE TIME =" $obs_time "s\n"
else
	obs_time=0
	echo -e "\tERROR: Couldn't get observation time from header, set obs_time "\
	        "to zero."
	exit
fi
if [ $(echo " $obs_time == 0.0" | bc) -eq 1 ]; then
	echo -e "\tERROR: Exposure time is zero. Exiting script."
	exit
fi

####################################
## Getting the mean energy spectrum
####################################

out_end="${out_root}_mean"
out_file="$out_dir/$out_end"

xspec_script="$out_dir/${prefix}_${day}_xspec.xcm"
spectrum_plot="${prefix}_${day}_meanonly"

if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"

cd "$out_dir"
if [ -e "${ccf_file}" ]; then
	python "$exe_dir"/energyspec.py "${ccf_file}" "${out_end}.${tab_ext}" \
			-b 0 -s 2
else
	echo -e "\tERROR: energyspec.py was not run. CCF output file does not "\
			"exist."
fi

if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ] && \
	[ -e "$bkgd_spec" ]; then

	cd "$out_dir"
	ascii2pha infile="${out_end}.${tab_ext}" \
		outfile="${out_end}.pha" \
		chanpres=yes \
		dtype=2 \
		qerror=yes \
		rows=- \
		tlmin=0 \
		detchans=64 \
		pois=no \
		telescope=XTE \
		instrume=PCA \
		detnam=PCU2 \
		filter=NONE \
		exposure=$obs_time \
		clobber=yes \
		respfile="$rsp_matrix" > $dump_file
# 	echo "XSPEC data: ${out_file}.pha"
else
	echo -e "\tERROR: ASCII2PHA was not run. Spectrum, response matrix, and/or"\
            " background spectrum do not exist."
fi
if [ ! -e "${out_end}.pha" ]; then
	echo -e "\tERROR: ASCII2PHA failed to create ${out_end}.pha."
fi
echo "data $out_end" >> $xspec_script
echo "ignore **-3.0 20.0-**" >> $xspec_script
echo "notice 3.0-20.0" >> $xspec_script
echo "ignore 11" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
echo "xsect vern" >> $xspec_script
echo "abund wilm" >> $xspec_script
echo "mod pow & 0 " >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@meanonly_points.pco 0.4 7 $spectrum_plot " >> $xspec_script
#echo "@meanonly_points.pco 0.2 0.7 $spectrum_plot " >> $xspec_script

echo "exit" >> $xspec_script

cd "$out_dir"
xspec < "$xspec_script" > $dump_file
if [ -e "$spectrum_plot.eps" ]; then open "$spectrum_plot.eps"; fi

################################################################################

#########################################################
## Making ccf deviation energy spectra and plotting them
#########################################################

spec_type=1  # 0 for mean+ccf, 1 for ccf, 2 for mean only
xspec_script="$out_dir/${prefix}_${day}_xspec.xcm"
spectrum_plot="${prefix}_${day}_ccfonly"

if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"
i=1

###############################################
## Generating energy spectra at each phase bin
###############################################

#for tbin in {6,13,19,24}; do
#for tbin in {20,25,30,35}; do
#for tbin in {8187,1,6,14}; do
for tbin in {8182,8187,1,6}; do

	out_end="${out_root}_ccfonly_${tbin}bin"

	## Make individual energy spectra from the CCF
	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py "${ccf_file}" \
		        "$out_dir/${out_end}.${tab_ext}" -b "$tbin" -s 1
	else
		echo -e "\tERROR: energyspec.py was not run. CCF output file does not" \
                "exist."
	fi

	## Put .dat spectra into .pha
	if [ -e "$rsp_matrix" ] && [ -e "$out_dir/${out_end}.${tab_ext}" ]; then

		cd "$out_dir"
		ascii2pha infile="${out_end}.${tab_ext}" \
			outfile="${out_end}.pha" \
			chanpres=yes \
			dtype=2 \
			qerror=yes \
			rows=- \
			tlmin=0 \
			detchans=64 \
			pois=no \
			telescope=XTE \
			instrume=PCA \
			detnam=PCU2 \
			filter=NONE \
			exposure=$obs_time \
			clobber=yes \
			respfile="$rsp_matrix" > $dump_file
	else
		echo -e "\tERROR: ASCII2PHA was not run. Spectrum and/or response" \
				"matrix do not exist."
	fi

	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tERROR: ASCII2pha did NOT run, so ${out_end}.pha does NOT "\
		        "exist."
	fi

	## Write to xspec script
	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

((i-=1))

#############################################
## Now we're ready to run xspec! -- ccf only
#############################################

echo "ignore 1-$i: **-3.0 20.0-**" >> $xspec_script
echo "notice 1-$i: 3.0-20.0" >> $xspec_script
echo "ignore 1-$i: 11" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
echo "xsect vern" >> $xspec_script
echo "abund wilm" >> $xspec_script
echo "mod pow+const*pow & 0 & & 0 & 0 & & " >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@ccfonly_points.pco -0.4 0.7 $spectrum_plot " >> $xspec_script
#echo "@ccfonly_points.pco -0.04 0.04 $spectrum_plot " >> $xspec_script
echo "exit" >> $xspec_script

cd "$out_dir"
#open "$xspec_script"
xspec < "$xspec_script" > $dump_file
if [ -e "$spectrum_plot.eps" ]; then open "$spectrum_plot.eps"; fi
cp "$spectrum_plot.eps" "$home_dir/Dropbox/Research/CCF_paper1/"

################################################################################

####################################################
## Making ccf+mean energy spectra and plotting them
####################################################

spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean
xspec_script="$out_dir/${prefix}_${day}_xspec.xcm"
spectrum_plot="${prefix}_${day}_ccfwmean"

if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"
i=1

###############################################
## Generating energy spectra at each phase bin
###############################################

# for tbin in {24,29,34,41}; do
#for tbin in {6,13,19,24}; do
#for tbin in {20,25,30,35}; do
#for tbin in {8187,1,6,14}; do
for tbin in {8182,8187,1,6}; do

	out_end="${out_root}_ccfwmean_${tbin}bin"

	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py "${ccf_file}" \
		        "$out_dir/${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\tERROR: energyspec.py was not run. CCF output file does not "\
                "exist."
	fi

	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ] && \
			[ -e "$bkgd_spec" ]; then

		cd "$out_dir"
		ascii2pha infile="${out_end}.${tab_ext}" \
			outfile="${out_end}.pha" \
			chanpres=yes \
			dtype=2 \
			qerror=yes \
			rows=- \
			tlmin=0 \
			detchans=64 \
			pois=no \
			telescope=XTE \
			instrume=PCA \
			detnam=PCU2 \
			filter=NONE \
			exposure=$obs_time \
			clobber=yes \
			respfile="$rsp_matrix" > $dump_file

	else
		echo -e "\tERROR: ASCII2PHA was not run. Spectrum, response matrix, "\
                "and/or background spectrum do not exist."
	fi

	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tERROR: ASCII2PHA failed to create ${out_end}.pha."
	fi

	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

((i-=1))
#############################################
## Now we're ready to run xspec! -- ccf+mean
#############################################

echo "ignore 1-$i: **-3.0 20.0-**" >> $xspec_script
echo "notice 1-$i: 3.0-20.0" >> $xspec_script
echo "ignore 1-$i: 11" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
echo "xsect vern" >> $xspec_script
echo "abund wilm" >> $xspec_script
echo "mod pow & 0 " >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@ccfwmean_points.pco 0.4 7 $spectrum_plot" >> $xspec_script
#echo "@ccfwmean_points.pco 0.2 0.7 $spectrum_plot" >> $xspec_script
echo "exit" >> $xspec_script

cd "$out_dir"
xspec < "$xspec_script" > $dump_file
if [ -e "$spectrum_plot.eps" ]; then open "$spectrum_plot.eps"; fi
cp "$spectrum_plot.eps" "$home_dir/Dropbox/Research/CCF_paper1/"

################################################################################