#!/bin/bash

home_dir=$(ls -d ~)  # the -d flag is extremely important here
in_dir="$home_dir/Dropbox/Research/cross_correlation"
exe_dir="$home_dir/Dropbox/Research/energy_spectra"
out_dir="$exe_dir/out_es"
script_dir="$home_dir/Dropbox/Scripts"
dump_file=dum.dat # Name of dumping file for intermediary steps

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi

propID=$1
obsID_list=$2
dt=$3
numsec=$4
testing=$5
day=$6

echo "$day"

# propID="P70080"
# obsID_list="$home_dir/Dropbox/Lists/${propID}_obsIDs.lst"
# dt=1
# numsec=4
# testing=0

obsID=$(head -1 $obsID_list)

# data_dir="$home_dir/Reduced_data/$propID"
data_dir="$home_dir/Dropbox/Research/sample_data"

ccf_file="$in_dir/out_ccf/${propID}_${day}_t${dt}_${numsec}sec.fits"
if (( testing==1 )); then
	ccf_file="$in_dir/out_ccf/test_${propID}_${day}_t${dt}_${numsec}sec.fits"
fi

tab_ext="dat"
plot_ext="png"

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. $HEADAS/headas-init.sh
fi

## Reading exposure time of observation
if [ -e "$ccf_file" ]; then
	if [ "${ccf_file##*.}" == "dat" ]; then
		obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$ccf_file')")
	elif [ "${ccf_file##*.}" == "fits" ]; then
		obs_time=$(python -c "from tools import get_key_val;  print get_key_val('$ccf_file', 0, 'EXPOSURE')")
	fi
	echo -e "EXPOSURE TIME =" $obs_time "s\n"
else
	obs_time=0
	echo -e "\tERROR: Couldn't get observation time from header, set obs_time to zero."
	exit
fi
if [ $(echo " $obs_time == 0.0" | bc) -eq 1 ]; then
	echo -e "\tERROR: Exposure time is zero. Exiting script."
	exit
fi

# rsp_dump_file="rsp_matrix_dump.dat"
# ## Making response matrix
# if [ -e "$rsp_matrix" ]; then
# 	echo "$rsp_matrix already exists."
# elif [ -e "$data_dir/all_evt.pha" ]; then
# 	pcarsp -f "$data_dir/all_evt.pha" -a "$data_dir/all.xfl" -l all -j y -p 2 -m n -n "$rsp_matrix" -z > $rsp_dump_file
# # 	pcarsp -f "$data_dir/all_evt.pha" -a "$data_dir/all.xfl" -l all -j y -p 2 -m n -n "$rsp_matrix" -z
# else
# 	echo -e "\tERROR: $data_dir/event.pha does NOT exist. pcarsp was NOT run."
# fi
# 
# # raw_event_bkgd="$data_dir/event_bkgd_notbinned.pha"
# raw_event_bkgd="$data_dir/evt_bkgd.pha"
# bkgd_spec="$out_dir/${propID}_${day}_event_spec.bkgd"
# 
# if [ -e "$raw_event_bkgd" ] && [ -e "$exe_dir/chan.txt" ] ; then
# 	rbnpha infile="$raw_event_bkgd" \
# 		outfile="$bkgd_spec" \
# 		binfile="$exe_dir/chan.txt" \
# 		clobber=yes
# else
# 	echo -e "\tERROR: $raw_event_bkgd and/or $exe_dir/chan.txt do NOT exist. rbnpha was NOT run."
# fi

bkgd_spec="$data_dir/evt_bkgd_rebinned.pha"
rsp_matrix="$data_dir/PCU2.rsp"
# rsp_matrix="$out_dir/${propID}_${day}_PCU2.rsp"
# quaternions="$data_dir/appx_quat.fits"
out_root="${propID}_${day}_t${dt}_${numsec}sec"


## Getting the mean energy spectrum
out_end="${out_root}_mean"
out_file="$out_dir/$out_end"

cd "$out_dir"
if [ -e "${ccf_file}" ]; then
	python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b 0 -s 2
else
	echo -e "\tERROR: ${ccf_file} does NOT exist, energyspec.py was NOT run."
fi
	
if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ] && [ -e "$bkgd_spec" ]; then
		
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
		respfile="$rsp_matrix" \
		backfile="$bkgd_spec" > $dump_file
# 	echo "XSPEC data: ${out_file}.pha"
else
	echo -e "\tSpectrum, response matrix, and/or background spectrum do NOT exist, ascii2pha was NOT run."
fi
if [ ! -e "${out_end}.pha" ]; then
	echo -e "\tERROR: ASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
fi


## Making ccf deviation spectra and plotting them
spec_type=1  # 0 for mean+ccf, 1 for ccf, 2 for mean only
xspec_script="$out_dir/${propID}_${day}_xspec.xcm"
spectrum_plot="${propID}_${day}_ccf"

if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi
touch "$xspec_script"
i=1

## Generating energy spectra at each phase bin
## for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
# for (( tbin=20; tbin<=35; tbin+=5 )); do
# # 	echo "$tbin"
# 	out_end="${out_root}_ccf_${tbin}bin"
# 
# 	if [ -e "${ccf_file}" ]; then
# 		python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
# 	else
# 		echo -e "\tERROR: ${ccf_file} does not exist, energyspec.py was NOT run."
# 	fi
# 	
# 	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
# 		
# 		cd "$out_dir"
# 		ascii2pha infile="${out_end}.${tab_ext}" \
# 			outfile="${out_end}.pha" \
# 			chanpres=yes \
# 			dtype=2 \
# 			qerror=yes \
# 			rows=- \
# 			tlmin=0 \
# 			detchans=64 \
# 			pois=no \
# 			telescope=XTE \
# 			instrume=PCA \
# 			detnam=PCU2 \
# 			filter=NONE \
# 			exposure=$obs_time \
# 			clobber=yes \
# 			respfile="$rsp_matrix" > $dump_file
# # 		echo "XSPEC data: ${out_file}.pha"
# 	else
# 		echo -e "\tSpectrum and/or response matrix do NOT exist, ascii2pha was NOT run."
# 	fi
# 	if [ ! -e "${out_end}.pha" ]; then
# 		echo -e "\tERROR: ASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
# 	fi
# 	
# 	echo "data $i:$i $out_end" >> $xspec_script
# 	((i+=1))
# done
# 
# ## Now we're ready to run xspec!
# echo "ignore 1-4: **-3 11 31-**" >> $xspec_script
# echo "notice 1-4: 3 31" >> $xspec_script
# echo "cpd /xw" >> $xspec_script
# echo "setplot energy" >> $xspec_script
# echo "systematic 0.005" >> $xspec_script
# echo "xsect vern" >> $xspec_script
# echo "abund wilm" >> $xspec_script
# echo "mod pow & 0" >> $xspec_script
# echo "iplot eeufspec" >> $xspec_script
# echo "@ccf_nomean.pco $spectrum_plot" >> $xspec_script
# echo "exit" >> $xspec_script
# cd "$out_dir"
# xspec < "$xspec_script" > "$dump_file"
# if [ -e "$spectrum_plot.eps" ]; then open "$spectrum_plot.eps"; fi


## Making ccf+mean spectra and plotting them
spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean
xspec_script="$out_dir/${propID}_${day}_xspec.xcm"
spectrum_plot="${propID}_${day}_ccfwmean"

if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi
touch "$xspec_script"
i=1

## Generating energy spectra at each phase bin
## for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=20; tbin<=35; tbin+=5 )); do
# 	echo "$tbin"
	out_end="${out_root}_ccfwmean_${tbin}bin"

	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\tERROR: ${ccf_file} does not exist, energyspec.py was NOT run."
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ] && [ -e "$bkgd_spec" ]; then
		
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
			respfile="$rsp_matrix" \
			backfile="$bkgd_spec" > $dump_file
# 		echo "XSPEC data: ${out_file}.pha"
	else
		echo -e "\tSpectrum, response matrix, and/or background spectrum do NOT exist, ascii2pha was NOT run."
	fi
	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tERROR: ASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
	fi
	
	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

## Now we're ready to run xspec!
echo "ignore 1-4: **-3 11 31-**" >> $xspec_script
echo "notice 1-4: 3 31" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
echo "xsect vern" >> $xspec_script
echo "abund wilm" >> $xspec_script
echo "mod pow & 0" >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@mean.pco $spectrum_plot" >> $xspec_script
echo "exit" >> $xspec_script
cd "$out_dir"
xspec < "$xspec_script" > "$dump_file"
if [ -e "$spectrum_plot.eps" ]; then open "$spectrum_plot.eps"; fi


###########################
## For testing out things
###########################


spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean
xspec_script="$out_dir/${propID}_${day}_xspec.xcm"
spectrum_plot="${propID}_${day}_ccftest"

if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi
touch "$xspec_script"
i=1

mod_vals=""
freeze_pars=""
j=1
k=2
l=4

## Generating energy spectra at each phase bin
for tbin in {20,27,33}; do  ## works for specific bin numbers
## for tbin in {20..35..5}; do  ## should work in bash 4.*, but i have bash 3.*
## for (( tbin=20; tbin<=35; tbin+=5 )); do  ## works for a range
# 	echo "$tbin"
	out_end="${out_root}_ccfwmean_${tbin}bin"

	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\tERROR: ${ccf_file} does not exist, energyspec.py was NOT run."
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ] && [ -e "$bkgd_spec" ]; then
		
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
			respfile="$rsp_matrix" \
			backfile="$bkgd_spec" > $dump_file
# 		echo "XSPEC data: ${out_file}.pha"
# 		echo -e "XSPEC resp: $rsp_matrix"
	else
		echo -e "\tSpectrum, response matrix, and/or background spectrum do NOT exist, ascii2pha was NOT run."
	fi
	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tERROR: ASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
	fi
	
	echo "data $i:$i $out_end" >> $xspec_script
# 	echo "data $i:$i $out_end"
	mod_vals+="& .113 & 1.9 .001 1.0 1.0 2.0 2.0 & .18 .001 .15 .15 .2 .2 & .8 .001 .6 .6 1. 1. & 50. .01 45. 45. 75. 75. "
	freeze_pars+="$j "
	((i+=1))
	((j+=5))
done
((i-=1))
# echo "$mod_vals"

## Now we're ready to run xspec!
echo "ignore 1-$i: **-3 11 31-**" >> $xspec_script
echo "notice 1-$i: 3 31" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
# echo "xsect vern" >> $xspec_script
# echo "abund wilm" >> $xspec_script
echo "mod phabs*(pow+bbodyrad) $mod_vals" >> $xspec_script
echo "freeze $freeze_pars" >> $xspec_script
echo "newpar 7 = 2" >> $xspec_script
echo "newpar 12 = 2" >> $xspec_script
echo "newpar 14 = 4" >> $xspec_script
echo "newpar 9 = 4" >> $xspec_script
echo "fit 1000" >> $xspec_script
echo "iplot resid" >> $xspec_script
open "$xspec_script"


## Now we're ready to run xspec!
# echo "ignore 1-$i: **-3 11 31-**" >> $xspec_script
# echo "notice 1-$i: 3 31" >> $xspec_script
# echo "cpd /xw" >> $xspec_script
# echo "setplot energy" >> $xspec_script
# echo "systematic 0.005" >> $xspec_script
# echo "xsect vern" >> $xspec_script
# echo "abund wilm" >> $xspec_script
# echo "mod pow & 0" >> $xspec_script
# echo "iplot eeufspec" >> $xspec_script
# echo "@ccf_points.pco $spectrum_plot" >> $xspec_script
# echo "exit" >> $xspec_script
# cd "$out_dir"
# xspec < "$xspec_script" > "$dump_file"
# if [ -e "$spectrum_plot.eps" ]; then open "$spectrum_plot.eps"; fi
