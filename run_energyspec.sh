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
day="141029"
dump_file=dum.dat # Name of dumping file for intermediary steps

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

propID=$1
obsID_list=$2
dt=$3
numsec=$4
testing=$5

# propID="P70080"
# obsID_list="$home_dir/Lists/${propID}_obsIDs.lst"
# dt=1
# numsec=4
# testing=0

obsID=$(head -1 $obsID_list)

ccf_file="$in_dir/out_ccf/${propID}_${day}_t${dt}_${numsec}sec.dat"
if (( testing==1 )); then
	ccf_file="$in_dir/out_ccf/test_${propID}_${day}_t${dt}_${numsec}sec.dat"
fi
data_dir="$home_dir/Reduced_data/$propID/$obsID"
# rsp_matrix="$data_dir/PCU2.rsp"
rsp_matrix="$out_dir/${propID}_${day}_${obsID}_PCU2.rsp"
quaternions="$data_dir/appx_quat.fits"

tab_ext="dat"
plot_ext="png"

## Making response matrix
if [ -e "$rsp_matrix" ]; then
	echo "$rsp_matrix already exists."
elif [ -e "$data_dir/event.pha" ] && [ -e "${quaternions}" ]; then
	pcarsp -f "$data_dir/event.pha" -a "${quaternions}" -l all -j y -p 2 -m n -n "$rsp_matrix"
else
	echo -e "\n\t $data_dir/event.pha and/or ${quaternions} do NOT exist. pcarsp was NOT run.\n"
fi

## Reading exposure time of observation
if [ -e "$ccf_file" ]; then
	if [ "${ccf_file##*.}" == "dat" ]; then
		obs_time=$(python -c "from tools import read_obs_time; print read_obs_time('$ccf_file')")
	elif [ "${ccf_file##*.}" == "fits" ]; then
		obs_time=$(python -c "from tools import get_key_val;  print get_key_val('$ccf_file', 0, 'EXPOSURE')")
	fi
	echo "EXPOSURE TIME =" $obs_time "s"
else
	obs_time=0
	echo -e "\tERROR: Couldn't get observation time from header, set obs_time to zero."
	exit
fi
if [ $(echo " $obs_time == 0.0" | bc) -eq 1 ]; then
	echo -e "\tERROR: Exposure time is zero. Exiting script."
	exit
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

out_root="${propID}_${day}_t${dt}_${numsec}sec"


## Getting the mean energy spectrum
out_end="${out_root}_mean"
out_file="$out_dir/$out_end"

cd "$out_dir"
if [ -e "${ccf_file}" ]; then
	python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b 0 -s 2
else
	echo -e "\t${ccf_file} does NOT exist, energyspec.py was NOT run."
fi
	
if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
		
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
		respfile="$rsp_matrix" > $dump_file
	echo "XSPEC data: ${out_file}.pha"
	echo -e "XSPEC resp: $rsp_matrix"
else
	echo -e "\t${rsp_matrix} and/or ${out_end}.${tab_ext} do NOT exist, ascii2pha was NOT run."
fi
if [ ! -e "${out_end}.pha" ]; then
	echo -e "\tASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
fi


## Now we're ready to run xspec!

## Making ccf+mean spectra and plotting them
spec_type=1  # 0 for mean+ccf, 1 for ccf, 2 for mean only
xspec_script="$out_dir/${propID}_${day}_xspec.xcm"
spectrum_plot="${propID}_${day}_ccf"

if [ -e "$xspec_script" ]; then
	rm "$xspec_script"
fi
touch "$xspec_script"
i=1

## Generating energy spectra at each phase bin
# for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=20; tbin<=35; tbin+=5 )); do
# 	echo "$tbin"
	out_end="${out_root}_ccf_${tbin}bin"

	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\t${ccf_file} does not exist, energyspec.py was NOT run."
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
		
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
			respfile="$rsp_matrix" > $dump_file
		echo "XSPEC data: ${out_file}.pha"
		echo -e "XSPEC resp: $rsp_matrix"
	else
		echo -e "\t${rsp_matrix} and/or ${out_end}.${tab_ext} do NOT exist, ascii2pha was NOT run."
	fi
	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
	fi
	
	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

echo "ignore 1-4: **-3 11 27-**" >> $xspec_script
echo "notice 1-4: 3 27" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "mod pow & 0" >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@ccf_nomean.pco $spectrum_plot" >> $xspec_script
echo "exit" >> $xspec_script
cd "$out_dir"
xspec < "$xspec_script" > "$dump_file"
if [ -e "$spectrum_plot.eps" ]; then
	open "$spectrum_plot.eps"
fi


## Making ccf+mean spectra and plotting them
spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean
xspec_script="$out_dir/${propID}_${day}_xspec.xcm"
spectrum_plot="${propID}_${day}_ccfwmean"

if [ -e "$xspec_script" ]; then
	rm "$xspec_script"
fi
touch "$xspec_script"
i=1

## Generating energy spectra at each phase bin
# for tbin in {25..40..5}; do  ## should work in bash 4.*, but i have 3.2.*
for (( tbin=20; tbin<=35; tbin+=5 )); do
# 	echo "$tbin"
	out_end="${out_root}_ccfwmean_${tbin}bin"

	if [ -e "${ccf_file}" ]; then
		python "$exe_dir"/energyspec.py -i "${ccf_file}" -o "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
	else
		echo -e "\t${ccf_file} does not exist, energyspec.py was NOT run."
	fi
	
	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
		
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
			respfile="$rsp_matrix" > $dump_file
		echo "XSPEC data: ${out_file}.pha"
		echo -e "XSPEC resp: $rsp_matrix"
	else
		echo -e "\t${rsp_matrix} and/or ${out_end}.${tab_ext} do NOT exist, ascii2pha was NOT run."
	fi
	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tASCII2pha did NOT run, so ${out_end}.pha does NOT exist."
	fi
	
	echo "data $i:$i $out_end" >> $xspec_script
	((i+=1))
done

echo "ignore 1-4: **-3 11 27-**" >> $xspec_script
echo "notice 1-4: 3 27" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "mod pow & 0" >> $xspec_script
echo "iplot eeufspec" >> $xspec_script
echo "@mean.pco $spectrum_plot" >> $xspec_script
echo "exit" >> $xspec_script
cd "$out_dir"
xspec < "$xspec_script" > "$dump_file"
if [ -e "$spectrum_plot.eps" ]; then
	open "$spectrum_plot.eps"
fi

## Now we're ready to run xspec!