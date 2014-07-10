#!/bin/bash

#####################################
## 
## Need to have heainit running!!
##
#####################################

home_dir=$(ls -d ~)  # the -d flag is extremely important here
in_dir="$home_dir/Dropbox/Research/cross_correlation"
exec_dir="$home_dir/Dropbox/Research/energy_spectra"
out_dir="$exec_dir/out_es"
script_dir="$home_dir/Dropbox/Scripts"
# day=$(date +%y%m%d)  # make the date a string and assign it to 'day', for filename
day="140624"

if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

propID=$1
obsID_list=$2
dt=$3
numsec=$4

obsID=$(head -1 $obsID_list)
filt_loc=1

in_file="$in_dir/out_ccf/${propID}_${day}_t${dt}_${numsec}sec.dat"
data_dir="$home_dir/Reduced_data/$propID/$obsID"
# rsp_matrix="$data_dir/PCU2.rsp"
rsp_matrix="$out_dir/${propID}_${day}_${obsID}_PCU2.rsp"
quaternions="$data_dir/appx_quat.fits"

tab_ext="dat"
plot_ext="png"

if [ -e "$rsp_matrix" ]; then
	echo "$rsp_matrix already exists."
elif [ -e "$data_dir/event.pha" ] && [ -e "${quaternions}" ]; then
	pcarsp -f "$data_dir/event.pha" -a "${quaternions}" -l all -j y -p 2 -m n -n "$rsp_matrix"
else
	echo -e "\n\t $data_dir/event.pha and/or ${quaternions} do NOT exist. pcarsp was NOT run.\n"
fi
if [ -e "$in_file" ]; then
	obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$in_file')")
	echo "OBS TIME = " $obs_time "seconds"
else
	obs_time = 0
	echo -e "\n\t Couldn't get observation time from header, set obs_time to zero.\n"
fi

for pbin in {25..26}; do

	out_end="${propID}_${day}_t${dt}_${numsec}sec_pbin_${pbin}"
	out_file="$out_dir/$out_end"

	if [ -e "$in_file" ]; then
		python "$exec_dir"/energyspec.py "$in_file" "${out_file}.${tab_ext}" "$pbin" "$filt_loc"
	else
		echo -e "\n\t $in_file does not exist, energyspec.py was NOT run.\n"
	fi
	
	pwd
	if [ -e "$rsp_matrix" ] && [ -e "${out_file}.${tab_ext}" ]; then
		
		cd "$out_dir"
		# ascii2pha infile="${out_end}.${tab_ext}" \
# 			outfile="${out_end}.pha" \
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
# 		echo "XSPEC data: ${out_file}.pha"
# 		echo -e "XSPEC resp: $rsp_matrix\n"
	else
		echo -e "\n\t ${rsp_matrix} and/or ${out_file}.${tab_ext} do NOT exist, ascii2pha was NOT run.\n"
	fi
	
	if [ -e "${out_file}.${tab_ext}" ]; then
		python "$exec_dir"/ratio_spectrum.py "${out_file}.${tab_ext}" "$obsID_list" "${out_file}_rwm.${tab_ext}"
	else
		echo -e "\n\t ${out_file}.${tab_ext} does not exist, ratio_spectrum.py was NOT run."
	fi
	if [ -e "$rsp_matrix" ] && [ -e "${out_file}_rwm.${tab_ext}" ]; then
		cd "$out_dir"
		ascii2pha infile="${out_end}_rwm.${tab_ext}" \
			outfile="${out_end}_rwm.pha" \
			chanpres=no \
			dtype=2 \
			qerror=yes \
			rows=- \
			fchan=0 \
			tlmin=0 \
			detchans=64 \
			telescope=RXTE \
			instrume=PCA \
			detnam=PCU2 \
			filter=NONE \
			exposure=$obs_time \
			clobber=yes \
			respfile="$rsp_matrix"
		echo "XSPEC data: ${out_file}_rwm.pha"
		echo -e "XSPEC resp: $rsp_matrix\n"
	else
		echo -e "\n\t ${rsp_matrix} and/or ${out_file}_rwm.${tab_ext} do NOT exist, ascii2pha was NOT run."
	fi

done

## Now we're ready to run xspec!