#! /bin/bash

#########################################################################################
##
## Make a mean energy spectrum of event-mode data for all obsIDs of a propID.
## 
## Notes: heainit needs to already be running!
## 
## 
#########################################################################################
day=$(date +%y%m%d) 
echo "Starting script 'mean_spectrum.sh', $day"
echo -e "If 'heainit' isn't running, you'll get errors!\n"

home_dir=$(ls -d ~)  # The home directory of this machine; the -d flag is extremely important here
script_dir="$home_dir/Dropbox/Scripts" # Directory containing the scripts
list_dir="$home_dir/Lists" # A folder of lists; tells which files we're using
propID=$1 # Proposal ID
obsID_list=$2 # List of relevant obsIDs for the given propID
out_dir="$home_dir/Reduced_data/$propID" # Output directory
dump_file=dum.dat # Name of dumping file for intermediary steps
all_evt="$list_dir/all_evt.lst" # List of all evt files, to make a mean spectrum
all_gti="$out_dir/all.gti"
all_filter="$out_dir/all.xfl"
filter_list="$list_dir/all_filter.lst"
o_filter_list="$list_dir/all_filter_ordered.lst"
gti_list="$list_dir/all_gti.lst"

if [ -e "$all_evt" ]; then
	rm "$all_evt"
fi
if [ -e "$filter_list" ]; then
	rm "$filter_list"
fi
if [ -e "$o_filter_list" ]; then
	rm "$o_filter_list"
fi
if [ -e "$gti_list" ]; then
	rm "$gti_list"
fi
touch "$all_evt"
touch "$filter_list"
touch "$o_filter_list"
touch "$gti_list"

## Make a list of all evt files and filter files
for obsID in $(cat $obsID_list); do
# 	echo "$obsID"
	ls $out_dir/$obsID/evt*.pca >> "$all_evt"
# 	echo "$out_dir/$obsID/filter.xfl"
	echo "$out_dir/$obsID/filter.xfl" >> "$filter_list"
	echo "$out_dir/$obsID/gti_file.gti" >> "$gti_list"
done

## Time-sort the filter files (switched quotes to avoid syntax error)
python -c "from tools import time_ordered_list; time_ordered_list('$filter_list')" >> "$o_filter_list"

## Merge filter files
fmerge infiles=@"$o_filter_list" outfile="$all_filter" columns=- clobber=yes
# echo "merged filter: $all_filter"
## Replace tstop, date-end and time-end
last_file=$(cat $o_filter_list | tail -n 1)
# echo "last file = $last_file"
last_tstop=$(python -c "from tools import get_key_val; print get_key_val('$last_file',0,'TSTOP')")
last_dateend=$(python -c "from tools import get_key_val; print get_key_val('$last_file',0,'DATE-END')")
last_timeend=$(python -c "from tools import get_key_val; print get_key_val('$last_file',0,'TIME-END')")

# echo "last tstop = $last_tstop"
# echo "last dateend = $last_dateend"
# echo "last timeend = $last_timeend"

python -c "from tools import replace_key_val; replace_key_val('$all_filter', 0, 'TSTOP', float('$last_tstop')); replace_key_val('$all_filter', 0, 'DATE-END', '$last_dateend'); replace_key_val('$all_filter', 0, 'TIME-END', '$last_timeend'); replace_key_val('$all_filter', 1, 'TSTOP', float('$last_tstop'))"


## Make GTI of mega-filter
bin_loc=$(python -c "from tools import get_key_val; print get_key_val('$all_filter', 0, 'TIMEPIXR')")
# echo "TIMEPIXR = $bin_loc"
filtex="(PCU2_ON==1)&&(PCU0_ON==1)&&(elv>10)&&(offset<0.02)&&(VpX1LCntPcu2<=150)&&(VpX1RCntPcu2<=150)" # same filter expression as used in gti_and_bkgd.sh
# not getting rid of pcu0 in here because i want the same times as the manual 
#  runs; bitfile is ensuring that it's only using pcu2

if (( bin_loc == 0 )); then
# 	echo "In here 1"
	maketime infile=$all_filter outfile=$all_gti expr=$filtex name=NAME \
		value=VALUE time=Time compact=no clobber=yes prefr=0.0 postfr=1.0
elif (( bin_loc == 1 )); then
# 	echo "In here 2"
	maketime infile=$all_filter outfile=$all_gti expr=$filtex name=NAME \
		value=VALUE time=Time compact=no clobber=yes prefr=1.0 postfr=0.0
else
# 	echo "In here 3"
	maketime infile=$all_filter outfile=$all_gti expr=$filtex name=NAME \
		value=VALUE time=Time compact=no clobber=yes prefr=0.5 postfr=0.5
fi

run_log="run.log"

#############################################################
## Making an event-mode spectrum for all files in the propID
#############################################################

# seextrct infile=@"$all_evt" \
# 	maxmiss=INDEF \
# 	gtiorfile=- \
# 	gtiandfile="$all_gti" \
# 	outroot="$out_dir/all_evt" \
# 	bitfile="$list_dir"/bitfile_evt_PCU2 \
# 	timecol="TIME" \
# 	columns="Event" \
# 	multiple=yes \
# 	binsz=1 \
# 	printmode=BOTH \
# 	lcmode=RATE \
# 	spmode=SUM \
# 	timemin=INDEF \
# 	timemax=INDEF \
# 	timeint=INDEF \
# 	chmin=INDEF \
# 	chmax=INDEF \
# 	chint=INDEF \
# 	chbin=INDEF \
# 	mode=ql > run.log
# 	
# if [ -e $dump_file ]; then
# 	rm -f $dump_file
# fi
# 
# if [ ! -e "$out_dir/all_evt.pha" ]; then
# 	echo -e "\n\tERROR: $out_dir/all_evt.pha not made!\n"
# fi
# 
# rsp_matrix="$out_dir/all_PCU2.rsp"
# first_obsID=$(cat $obsID_list | head -n 1)
# # echo "$first_obsID"
# quaternions="$out_dir/$first_obsID/appx_quat.fits"
# 
# 
# if [ -e "$rsp_matrix" ]; then
# 	echo "$rsp_matrix already exists."
# elif [ -e "$out_dir/all_evt.pha" ] && [ -e "${quaternions}" ]; then
# 	pcarsp -f "$out_dir/all_evt.pha" -a "${quaternions}" -l all -j y -p 2 -m n -n "$rsp_matrix"
# else
# 	echo -e "\n\t $data_dir/event.pha and/or ${quaternions} do NOT exist. pcarsp was NOT run.\n"
# fi
# 
# # python "$script_dir/plot_evt_lc.py" "$propID" "$out_dir/all_evt.lc" lc_plot.png
# 
# seextrct_exposure=$(python -c "from tools import get_key_val; print get_key_val('$out_dir/all_evt.pha',1,'EXPOSURE')")
# echo "SEEXTRCT EXPOSURE TIME: $seextrct_exposure s"
# echo `echo "$seextrct_exposure / 16.0" | bc -l`