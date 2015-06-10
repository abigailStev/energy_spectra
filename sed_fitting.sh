#!/bin/bash

################################################################################
## 
## Bash script for phase-resolved spectroscopy: run energyspec.py to make phase-
## resolved energy spectra, make an XSPEC energy spectrum fitting script, run 
## the script, read off fit data from log file with multispec_plots.py, and make
## plots of fit parameters changing with QPO phase.
##
## Example call: ./sed_fitting.sh GX339-BQPO 64 64 0 150131
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries) 
## 		  must be installed in order to run this script. 
##
## Written by Abigail Stevens, A.L.Stevens at uva.nl, 2015
##
################################################################################

## Checking the number of input arguments
if (( $# != 5 )); then
    echo -e "\tUsage: ./sed_fitting.sh <prefix> <dt multiple> <num "\
        "seconds> <testing> <date>\n"
    exit
fi

## ./sed_fitting.sh GX339-BQPO 64 64 0 150526
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
out_dir="$exe_dir/out_es"
dump_file=dump.txt # Name of dumping file for intermediary steps

data_dir="$home_dir/Reduced_data/$prefix"
# data_dir="$home_dir/Dropbox/Research/sample_data"

bkgd_spec="$data_dir/evt_bkgd_rebinned.pha"
rsp_matrix="$data_dir/PCU2.rsp"
out_root="${prefix}_${day}_t${dt}_${numsec}sec"

spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean

xspec_script="$out_dir/${prefix}_${day}_xspec.xcm"
spectrum_plot="${prefix}_${day}_allspectra"

ccf_file="$in_dir/out_ccf/${prefix}_150401_t${dt}_${numsec}sec.fits"
if (( $testing==1 )); then
	ccf_file="$in_dir/out_ccf/test_${prefix}_150401_t${dt}_${numsec}sec.fits"
fi

tab_ext="dat"
plot_ext="png"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"

obs_time=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'EXPOSURE')")

i=1
mod_vals=""

cd "$out_dir"

###############################################
## Generating energy spectra at each phase bin
###############################################

#for tbin in {6,13,19,24}; do
 for (( tbin=6; tbin<=30; tbin++ )); do

	out_end="${out_root}_ccfwmean_${tbin}bin"

#	if [ -e "${ccf_file}" ]; then
#		python "$exe_dir"/energyspec.py "${ccf_file}" "${out_end}.${tab_ext}" \
#			-b "$tbin" -s "$spec_type"
#	else
#		echo -e "\tERROR: ${ccf_file} does not exist, energyspec.py was NOT run."
#	fi
#
#	if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ] && \
#		[ -e "$bkgd_spec" ]; then
#		ascii2pha infile="${out_end}.${tab_ext}" \
#			outfile="${out_end}.pha" \
#			chanpres=yes \
#			dtype=2 \
#			qerror=yes \
#			rows=- \
#			tlmin=0 \
#			detchans=64 \
#			pois=no \
#			telescope=XTE \
#			instrume=PCA \
#			detnam=PCU2 \
#			filter=NONE \
#			exposure=$obs_time \
#			clobber=yes \
#			respfile="$rsp_matrix" \
#			backfile="$bkgd_spec" > $dump_file
#	else
#		echo -e "\tERROR: ASCII2PHA was not run. Spectrum, response matrix, "\
#            "and/or background spectrum do not exist."
#	fi

	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tERROR: ASCII2PHA failed to create ${out_end}.pha."
		echo -e "\tExiting script."
		exit
	fi
	
	## Writing file names to xspec script
	echo "data $i:$i $out_end.pha" >> $xspec_script
	
	## for phabs*(simpl*diskbb+gauss)
# 	mod_vals+="& & & 0.2 & & & & & & "  ## FracSctr
# 	mod_vals+="& & 3.0 & 0.2 & & & & & & " ## Gamma and FracSctr
#	mod_vals+="& & & 0.2 & & 0.8 & & & & " ## FracSctr and Tin
# 	mod_vals+="& & & 0.2 & & & 3000 & & & " ## FracSctr and norm(BB)
# 	mod_vals+="& & & 0.2 & & & & 6.33 .02 6.1 6.1 6.9 6.9 & & " ## FracSctr and lineE
# 	mod_vals+="& & & 0.2 & & & & & 0.8 & " ## FracSctr and Sigma
# 	mod_vals+="& & & 0.2 & & & & & & 0.01"  ## FracSctr and norm(E)
# 	mod_vals+="& & 3.0 & 0.2 & & 0.8 & & & & "  ## Gamma, FracSctr, and Tin
# 	xspec_log="${prefix}_${day}_xspecfit_G-FS-T.log"
# 	mod_vals+="& & 3.0 & 0.2 & & & 3000 & & & "  ## Gamma, FracSctr, and norm(BB)
#	mod_vals+="& & 3.0 & 0.2 & & 0.8 & 3000 & & & "  ## Gamma, FracSctr, Tin, and norm(BB)

## For phabs*(simpl*const*diskbb+diskbb+bbodyrad+gauss)
#	mod_vals+=" &  & &  &  & &   &  &  &  &  & & & & "
#	mod_vals+=" &  & 2.8 &  &  & .2 &   &  &  &  &  & & & & "
#	mod_vals+=" &  & 2.8 &  &  & .2 &  &  &  &  & .6 & & & & "
#	mod_vals+=" &  & 2.8 &  &  &    &   & &  &  &  & & & & "
#	mod_vals+=" &  & &  &  & &   &  & .8 &  &  & & & & "
#	mod_vals+=" &  & &  &  & &   &  &  & 3000 &  & & & & "
#	mod_vals+=" &  & &  &  & &   &  &  &  & .8 & & & & "
#	mod_vals+=" &  & &  &  & &   &  &  &  &  & 3000 & & & "


## For phabs*(simpl*bbodyrad+bbodyrad+gauss)
#	mod_vals+=" &  &  &     &  & &  &  &  &  & &  "
#	mod_vals+=" &  &  & 0.2 &  & &  &  &  &  & &   "
#	mod_vals+=" &  & 2.8 0.01 2.5 2.5 3.2 3.2 & &  & &  &  &  &  & &  "
#	mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 & 0.2 &  & & &  &  &  & &  "
#	mod_vals+=" &  & 2.8 0.01 2.5 2.5 3.2 3.2 &  &  & &  & .6 .002 0.1 0.1 0.9 0.9 &  &  & &  "
#	mod_vals+=" &  & 2.8 0.01 2.5 2.5 3.2 3.2 &  &  & &  &  & 2500 &  & &  "
#	mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 & 0.2 &  & & & .6 .002 0.1 0.1 0.9 0.9 &  &  & &  "

## For phabs*(nthcomp+bbodyrad+bbodyrad+gauss)
#	mod_vals+=" &  &     &  &  & &  & &  &  &  & & & & "
#	mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 &  &  & &  & &  &  &  & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 &  &  &  & & & & "
#    mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 &  &  & &  & .2 &  &  &  & & & & "
#	mod_vals+=" &  &  2.8 0.01 2.0 2.0 3.1 3.1  &  &  & &  & .2 &  &  & .6 .005 0.1 0.1 2.0 2.0 & & & & "
#	mod_vals+=" &  &  2.8 0.01 2.6 2.6 3.1 3.1  &  &  & &  & &  &  & .6 .005 0.1 0.1 2.0 2.0 & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 & .6 .005 0.1 0.1 2.0 2.0 &  &  & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 &    & 3000 &  & & & & "

	xspec_log="${prefix}_${day}_xspecfit_NTH-2BB-regerr.log"
	
	((i+=1))
done
((i-=1))

numpar=14
################################################################################

if [ -e "$xspec_log" ]; then rm "$xspec_log"; fi; touch "$xspec_log"
if [ -e "MCMC_NTH-2BB.fits" ]; then rm "MCMC_NTH-2BB.fits"; fi;

########################################
## Writing the rest of the xspec script
########################################

echo "ignore 1-$i: **-3.0 20.0-**" >> $xspec_script
echo "notice 1-$i: 3.0-20.0" >> $xspec_script
echo "ignore 1-$i: 11" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
echo "xsect vern" >> $xspec_script
echo "abund wilm" >> $xspec_script


#echo "mod phabs*(simpl*diskbb+bbodyrad+gauss) $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.8 0.01 2.6 2.6 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1" >> $xspec_script
#echo "freeze 4" >> $xspec_script
#echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
##echo "newpar 5 0.744958" >> $xspec_script
##echo "newpar 5 0.80" >> $xspec_script
##echo "freeze 5" >> $xspec_script
##echo "newpar 6 3214.06" >> $xspec_script
#echo "newpar 6 3300" >> $xspec_script
##echo "freeze 6" >> $xspec_script
#echo "newpar 7 0.6 0.002 0.1 0.1 0.9 0.9" >> $xspec_script
#echo "newpar 8 2500" >> $xspec_script
#echo "newpar 9 6.33 0.005 6.1 6.1 6.9 6.9" >> $xspec_script
#echo "newpar 10 0.8 .005 0.65 0.65 0.9 0.9" >> $xspec_script
#echo "newpar 11 1.0E-02" >> $xspec_script


echo "mod phabs*(nthcomp+bbodyrad+bbodyrad+gauss) $mod_vals" >> $xspec_script
echo "newpar 1 0.6" >> $xspec_script  ## This value is from Reynolds and Miller 2013
echo "freeze 1" >> $xspec_script
echo "newpar 2 2.8 0.01 2.0 2.0 3.1 3.1" >> $xspec_script
echo "newpar 3 80.0 50.0 50.0 200.0 200.0" >> $xspec_script
echo "thaw 3" >> $xspec_script
echo "newpar 4 0.8" >> $xspec_script
echo "thaw 4" >> $xspec_script
echo "newpar 5 0" >> $xspec_script
echo "freeze 5" >> $xspec_script
echo "newpar 6 0.0" >> $xspec_script
echo "freeze 6" >> $xspec_script
echo "newpar 7 0.5" >> $xspec_script
echo "newpar 8 .8 .005 0.5 0.5 1.0 1.0" >> $xspec_script
echo "newpar 9 3300 " >> $xspec_script
#echo "freeze 9" >> $xspec_script
echo "newpar 10 .6 .005 0.1 0.1 2.0 2.0" >> $xspec_script
echo "newpar 11 1.14574E+04" >> $xspec_script
echo "freeze 11" >> $xspec_script
#echo "newpar 11 3340.2592 " >> $xspec_script
echo "newpar 12 6.33 0.02 6.1 6.1 6.9 6.9" >> $xspec_script
echo "newpar 13 0.8 .005 0.65 0.65 1.0 1.0" >> $xspec_script
echo "newpar 14 1.0E-02" >> $xspec_script
## Tieing disk-component blackbodies together
for (( x=1; x<=i*numpar; x++ )); do
    n=$((x%numpar))
    y=$((x+4))
    if (( n == 4 )) && (( n != 0 )); then
        echo "newpar $x = $y" >> $xspec_script
    fi
done

#echo "mod phabs*(simpl*const*bbodyrad+bbodyrad+bbodyrad+gauss) $mod_vals" >> $xspec_script
# echo "mod phabs*(simpl*diskbb+gauss) $mod_vals" >> $xspec_script

#echo "newpar 1 0.6" >> $xspec_script  ## This value is from Reynolds and Miller 2013
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.8" >> $xspec_script
##echo "newpar 3 80.0 50.0 50.0 130.0 130.0" >> $xspec_script
##echo "thaw 3" >> $xspec_script
##echo "newpar 3 100.0" >> $xspec_script
##echo "freeze 3" >> $xspec_script
#echo "newpar 3 1" >> $xspec_script
#echo "freeze 3" >> $xspec_script
##echo "newpar 4 0.8" >> $xspec_script
##echo "thaw 4" >> $xspec_script
#echo "newpar 4 1" >> $xspec_script
#echo "freeze 4" >> $xspec_script
##echo "newpar 5 0" >> $xspec_script
#echo "newpar 5 0.2" >> $xspec_script
##echo "newpar 6 0.0" >> $xspec_script
#echo "newpar 6 0.8 0.005 0.5 0.5 1.0 1.0" < $xspec_script
##echo "newpar 7 1.0" >> $xspec_script
#echo "newpar 7 3300 " >> $xspec_script
#echo "freeze 7" >> $xspec_script
##echo "newpar 8 0.776372" >> $xspec_script
##echo "freeze 8" >> $xspec_script
#echo "newpar 8 .8 .005 0.5 0.5 1.0 1.0" >> $xspec_script
#echo "newpar 9 3300 " >> $xspec_script
#echo "freeze 9" >> $xspec_script
#echo "newpar 10 .6 .005 0.1 0.1 2.0 2.0" >> $xspec_script
#echo "newpar 11 2500" >> $xspec_script
#echo "freeze 11" >> $xspec_script
##echo "newpar 11 3340.2592 " >> $xspec_script
#echo "newpar 12 6.33 0.02 6.1 6.1 6.9 6.9" >> $xspec_script
#echo "newpar 13 0.791920" >> $xspec_script
#echo "newpar 14 1.68325E-02" >> $xspec_script

## Tieing disk-component blackbodies together
#for (( x=1; x<=i*numpar; x++ )); do
#    n=$((x%numpar))
#    y=$((x-2))
#    if (( n == 8 )) && (( n != 0 )); then
#        echo "newpar $x = $y" >> $xspec_script
#    elif (( n == 9 )) && (( n != 0 )); then
#        echo "newpar $x = $y" >> $xspec_script
#    fi
#done

echo "log $xspec_log" >> $xspec_script
#echo "chatter 4" >> $xspec_script
echo "fit 1000" >> $xspec_script
#echo "chain type gw " >> $xspec_script
#echo "chain burn 2000" >> $xspec_script
#echo "chain walkers 1000" >> $xspec_script
#echo "chain length 100000" >> $xspec_script
#echo "chain run MCMC_NTH-2BB.fits">> $xspec_script
#echo "error maximum 10000. 2.706 1-350" >> $xspec_script
echo "cpd /xw" >> $xspec_script
echo "iplo ratio" >> $xspec_script
echo "hardcopy ratio.eps/cps" >> $xspec_script
echo "quit" >> $xspec_script
# echo "show par" >> $xspec_scriptyy
cd "$out_dir"
xspec - "$xspec_script"
# > $dump_file
#open ratio.eps
#open "$xspec_log"
echo "Finished running sed_fitting.sh"
cd ..
pwd
#echo "$xspec_log"

#python multifit_plots.py "$out_dir/$xspec_log"

################################################################################