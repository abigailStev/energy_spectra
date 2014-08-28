#!/bin/bash

#####################################
## 
## Need to have heainit running!!
##
#####################################

home_dir=$(ls -d ~)  # the -d flag is extremely important here
in_dir="$home_dir/Dropbox/Research/cross_correlation"
exe_dir="$home_dir/Dropbox/Research/energy_spectra"
out_dir="$exe_dir/out_es"
script_dir="$home_dir/Dropbox/Scripts"
# day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
day="140826"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

propID=$1
obsID_list=$2
dt=$3
numsec=$4

obsID=$(head -1 $obsID_list)

ccf_file="$in_dir/out_ccf/${propID}_${day}_t${dt}_${numsec}sec.dat"
data_dir="$home_dir/Reduced_data/$propID/$obsID"
# rsp_matrix="$data_dir/PCU2.rsp"
rsp_matrix="$out_dir/${propID}_${day}_${obsID}_PCU2.rsp"
quaternions="$data_dir/appx_quat.fits"

tab_ext="dat"
plot_ext="png"

## Making response matrix, getting exposure time of observation
if [ -e "$rsp_matrix" ]; then
	echo "$rsp_matrix already exists."
elif [ -e "$data_dir/event.pha" ] && [ -e "${quaternions}" ]; then
	pcarsp -f "$data_dir/event.pha" -a "${quaternions}" -l all -j y -p 2 -m n -n "$rsp_matrix"
else
	echo -e "\n\t $data_dir/event.pha and/or ${quaternions} do NOT exist. pcarsp was NOT run.\n"
fi
if [ -e "$ccf_file" ]; then
	obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$ccf_file')")
	echo "OBS TIME =" $obs_time "seconds"
else
	obs_time = 0
	echo -e "\n\t Couldn't get observation time from header, set obs_time to zero.\n"
fi

## Making file with list of event.pha spectra
# spectra_list="$out_dir/${propID}_${day}_t${dt}_${numsec}sec_evtspectra.lst"
# if [ -e "$spectra_list" ]; then
# 	rm "$spectra_list"
# fi
# touch "$spectra_list"
# echo "$spectra_list"
# 
# for obsID in $(cat $obsID_list); do
# 	spec_file="$home_dir/Reduced_data/$propID/$obsID/event.pha"
# # 	echo "$spec_file"
# 	if [ -e "$spec_file" ]; then
#  		echo "$spec_file" >> $spectra_list
#  	fi
# done
# ## Mean spectrum
# ccf_plus_mean="$out_dir/${propID}_${day}_t${dt}_${numsec}sec_ccfspec"
# python "$exe_dir"/mean_spectrum.py -i "$spectra_list" -o "${ccf_plus_mean}.${tab_ext}" -c "$ccf_file"


## Generating energy spectra at each phase bin
# for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=25; tbin<=40; tbin+=5 )); do
# 	echo "$tbin"
	out_end="${propID}_${day}_t${dt}_${numsec}sec_pbin_${tbin}"
	out_file="$out_dir/$out_end"

	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_file}.${tab_ext}" -b "$tbin"
	else
		echo -e "\n\t ${ccf_file}.${tab_ext} does not exist, energyspec.py was NOT run.\n"
	fi
	
	pwd
	if [ -e "$rsp_matrix" ] && [ -e "${out_file}.${tab_ext}" ]; then
		
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
			telescope=RXTE \
			instrume=PCA \
			detnam=PCU2 \
			filter=NONE \
			exposure=$obs_time \
			clobber=yes \
			respfile="$rsp_matrix"
		echo "XSPEC data: ${out_file}.pha"
		echo -e "XSPEC resp: $rsp_matrix\n"
	else
		echo -e "\n\t${rsp_matrix} and/or ${out_file}.${tab_ext} do NOT exist, ascii2pha was NOT run.\n"
	fi

## For running ratio_spectrum

# 	if [ -e "${out_file}.${tab_ext}" ]; then
# 		python "$exe_dir"/ratio_spectrum.py "${out_file}.${tab_ext}" "$obsID_list" "${out_file}_rwm.${tab_ext}"
# 	else
# 		echo -e "\n\t ${out_file}.${tab_ext} does not exist, ratio_spectrum.py was NOT run."
# 	fi
# 	if [ -e "$rsp_matrix" ] && [ -e "${out_file}_rwm.${tab_ext}" ]; then
# 		cd "$out_dir"
# 		ascii2pha infile="${out_end}_rwm.${tab_ext}" \
# 			outfile="${out_end}_rwm.pha" \
# 			chanpres=no \
# 			dtype=2 \
# 			qerror=yes \
# 			rows=- \
# 			fchan=0 \
# 			tlmin=0 \
# 			detchans=64 \
# 			telescope=RXTE \
# 			instrume=PCA \
# 			detnam=PCU2 \
# 			filter=NONE \
# 			exposure=$obs_time \
# 			clobber=yes \
# 			respfile="$rsp_matrix"
# 		echo "XSPEC data: ${out_file}_rwm.pha"
# 		echo -e "XSPEC resp: $rsp_matrix\n"
# 	else
# 		echo -e "\n\t ${rsp_matrix} and/or ${out_file}_rwm.${tab_ext} do NOT exist, ascii2pha was NOT run."
# 	fi

done

## Now we're ready to run xspec!