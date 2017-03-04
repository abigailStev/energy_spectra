#!/bin/bash

################################################################################
## 
## Bash script for phase-resolved spectroscopy: run energyspec.py to make phase-
## resolved energy spectra, make an XSPEC spectral fitting script, run
## the script, read off fit data from log file with multispec_plots.py, and make
## plots of spectral parameters changing with QPO phase, fit a function to those
## parameter variations.
##
## Example call: source ./sed_fitting.sh GX339-BQPO 64 64 0 150131
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries) 
## 		  must be installed in order to run this script. 
##
## Author: Abigail Stevens <A.L.Stevens at uva.nl> 2015-2016
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
out_dir="$exe_dir/out_es/${prefix}"
dump_file=dump.txt # Name of dumping file for intermediary steps

data_dir="$home_dir/Reduced_data/$prefix"
# data_dir="$home_dir/Dropbox/Research/sample_data"

bkgd_spec="$data_dir/evt_bkgd_rebinned.pha"
rsp_matrix="${prefix}_PCU2.rsp"
#out_name="${prefix}_${day}_t${dt}_${numsec}sec"
out_name="${prefix}_${day}_t${dt}_${numsec}sec_adj"

mcmc_flag=0  # 0 for no, 1 for yes

#fit_specifier="3pd-1BB-FS-G-Tin"
fit_specifier="pBB-FS-G-p"
#fit_specifier="1BB-FS-G-NE"
#fit_specifier="2BB-FS-G-kT"
#fit_specifier+="-fzs-fzNbb8857-2"
fit_specifier+="-fzs"
fit_specifier+="-fzNbb"
#fit_specifier+="-fzTin"

if (( mcmc_flag == 1 )); then
    fit_specifier+="-wMCMC"
fi

export fit_specifier

xspec_script="$out_dir/${prefix}_${day}_${fit_specifier}_xspec.xcm"
xspec_log="${prefix}_${day}_${fit_specifier}_xspec.log"
ratio_plot="${prefix}_${day}_${fit_specifier}_ratio"
first_plot="${prefix}_${day}_${fit_specifier}_firstspec"
spectrum_plot="${prefix}_${day}_${fit_specifier}_allspectra"
tex_tab_file="$home_dir/Dropbox/Research/CCF_paper1/spec_models_testing.txt"
parfit_file="$out_dir/${prefix}_${day}_${fit_specifier}_funcfit.txt"

ccf_file="$in_dir/out_ccf/${prefix}/${out_name}.fits"
if (( $testing==1 )); then
	ccf_file="$in_dir/out_ccf/${prefix}/test_${out_name}.fits"
fi

tab_ext="dat"
plot_ext="eps"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
#if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"
if [ ! -e "$parfit_file" ]; then touch "$parfit_file"; fi

obs_time=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 1, 'EXPOSURE')")
#echo "$obs_time"
echo "lmod simpler /Users/abigailstevens/Dropbox/Research/xspecmods/" >> $xspec_script

i=1
mod_vals=""

cd "$out_dir"

###############################################################
## Generating a spectral energy distribution at each phase bin
###############################################################

#for tbin in {6,13,19,24}; do
#for (( tbin=6; tbin<=30; tbin++ )); do

n_spectra=$(( 8205-8182+1 ))
#n_spectra=71
#echo "$n_spectra"
export n_spectra

for (( timebin=8182; timebin<=8205; timebin++ )); do
#for (( timebin=8159; timebin<=8229; timebin++ )); do

    if (( timebin>=8192 )); then
        tbin=$((timebin-8192))
    else
        tbin=$((timebin))
    fi

	out_end="${out_name}_ccfwmean_${tbin}bin"
    if [ ! -e "${out_end}.pha" ]; then
        if [ -e "${ccf_file}" ]; then
            python "$exe_dir"/energyspec.py "${ccf_file}" \
                    "${out_end}.${tab_ext}" -b "$tbin" -s 0
        else
            echo -e "\tERROR: ${ccf_file} does not exist, energyspec.py was NOT run."
        fi

        if [ -e "$rsp_matrix" ] && [ -e "${out_end}.${tab_ext}" ]; then
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
## Don't need to include background file since I've subtracted background
## count rate per energy channel from the mean count rate per CI
        else
            echo -e "\tERROR: ASCII2PHA was not run. Spectrum, response matrix"\
                    ", and/or background spectrum do not exist."
        fi
    fi

	if [ ! -e "${out_end}.pha" ]; then
		echo -e "\tERROR: ASCII2PHA failed to create ${out_end}.pha."
		echo -e "\tExiting script."
		exit
	fi
	
	## Writing file names to xspec script
	echo "data $i:$i $out_end.pha" >> $xspec_script

    ###########################################################################
    ## Go through and uncomment the two lines for the untied parameter
    ## combination you want. First line is what's fed into XSPEC ("mod_vals"),
    ## second line is what's written in the table ("varpar").
    ###########################################################################

##
## For phabs*(simpler*diskbb+gauss)
##
#	mod_vals+="& & & & & & & & & "  ##  all tied
#	varpar=" - "
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & & & & & & & " ## Gamma
#	varpar="\$\\Gamma\$"
#	mod_vals+="& & & 0.2 & & & & & & "  ## Fscatt
#	varpar="\$F_{\\text{scatt}}\$"
#	mod_vals+="& & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & " ## Tdisk
#	varpar="\$T_{\\text{disk}}\$"
#	mod_vals+="& & & & & & 3000 & & & " ##  norm(BB)
#	varpar="\$N_{\\text{disk}}\$"
#	mod_vals+="& & & & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## Eline
#	varpar="\$E_{\\text{line}}\$"
#	mod_vals+="& & & & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & " ## sigma
#	varpar="\$\\sigma_{\\text{line}}\$"
#	mod_vals+="& & & & & & & & & 0.01"  ## Nline
#	varpar="\$N_{\\text{line}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & & & & " ## Fscatt and Gamma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$"
#	mod_vals+="& & & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & " ## Fscatt and Tdisk
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{disk}}\$"
#	mod_vals+="& & & 0.2 & & & 3000 & & & " ## Fscatt and norm(BB)
#	varpar="\$F_{\\text{scatt}}\$, \$N_{\\text{disk}}\$"
#	mod_vals+="& & & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## Fscatt and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$E_{\\text{line}}\$"
#	mod_vals+="& & & 0.2 & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & " ## Fscatt and sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\sigma_{\\text{line}}\$"
#	mod_vals+="& & & 0.2 & & & & & & 0.01"  ## Fscatt and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$N_{\\text{line}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & "  ## Fscatt, Gamma, Tdisk
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{disk}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & 3000 & & & "  ## Fscatt, Gamma, Ndisk
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{disk}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## Fscatt, Gamma, Eline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$E_{\\text{line}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & "  ## Fscatt, Gamma, sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$\\sigma_{\\text{line}}\$"
# 	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & & & & 0.01"  ## Fscatt, Gamma, Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$"
#    mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & 3000 & & & "  ## Fscatt, Gamma, Tdisk, Ndisk
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{disk}}\$, \$N_{\\text{disk}}\$"
#    mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & 6.4 0.005 5.5 5.5 7.0 7.0 & & "  ## Fscatt, Gamma, Tdisk, Eline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{disk}}\$, \$E_{\\text{line}}\$"
#    mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & 0.8 0.002 0.3 0.3 1.0 1.0 & "  ## Fscatt, Gamma, Tdisk, sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{disk}}\$, \$\\sigma_{\\text{line}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & 0.01 "  ## Gamma, Fscatt, Tdisk, Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{disk}}\$, \$N_{\\text{line}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & 3000 & & & 0.1 "  ## Gamma, Fscatt, Nline, and Ndisk
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$, \$N_{\\text{disk}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & 0.01 " ## Fscatt, Gamma, Nline, Eline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$, \$E_{\\text{line}}\$"
#	mod_vals+="& & 2.6 0.005 1.0 1.0 3.1 3.1 & 0.2 & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & 0.01 " ## Fscatt, Gamma, Nline, sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$, \$\\sigma_{\\text{line}}\$"

##
## For phabs*(simpl*const*diskbb+diskbb+bbodyrad+gauss)
##
#	mod_vals+=" &  & &  &  & &   &  &  &  &  & & & & "
#	mod_vals+=" &  & 2.6 &  &  & .2 &   &  &  &  &  & & & & "
#	mod_vals+=" &  & 2.6 &  &  & .2 &  &  &  &  & .6 & & & & "
#	mod_vals+=" &  & 2.6 &  &  &    &   & &  &  &  & & & & "
#	mod_vals+=" &  & &  &  & &   &  & .8 &  &  & & & & "
#	mod_vals+=" &  & &  &  & &   &  &  & 3000 &  & & & & "
#	mod_vals+=" &  & &  &  & &   &  &  &  & .8 & & & & "
#	mod_vals+=" &  & &  &  & &   &  &  &  &  & 3000 & & & "

##
## For phabs*(simpler*diskbb+bbodyrad+gauss) or phabs*(simpler*bbodyrad+bbodyrad+gauss)
##
#	mod_vals+=" &  &  &     &  & &  &  &  &  & &  "  ## All tied
#	varpar=" - "
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & &  & &  &  &  &  & &  "  ## Gamma
#	varpar="\$\\Gamma\$"
#	mod_vals+=" &  &  & 0.2 &  & &  &  &  &  & &   "   ## Fscatt
#	varpar="\$F_{\\text{scatt}}\$"
#	mod_vals+=" &  &  &  &  & &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## Tbb
#	varpar="\$T_{\\text{bb}}\$"
#	mod_vals+=" &  & &  &  & &  &  & 400 &  & &  "  ## Nbb
#	varpar="\$N_{\\text{bb}}\$"
#	mod_vals+=" &  & &  &  & &  &  &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Eline
#	varpar="\$E_{\\text{line}}\$"
#	mod_vals+=" &  & &  &  & &  &  &  &  & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## sigma
#	varpar="\$\\sigma_{\\text{line}}\$"
#    mod_vals+=" &  & &  &  & &  &  &  & &  & 0.01 "  ## Nline
#	varpar="\$N_{\\text{line}}\$"
#    mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & &  &  &  & &  "  ## Fscatt and Gamma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$"
#	mod_vals+=" &  &  & 0.2 &  & &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## Fscatt and Tbb
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{bb}}\$"
#	mod_vals+=" &  &  & 0.2 &  & &  &  & 2500 &  & &  "  ## Fscatt and Nbb
#	varpar="\$F_{\\text{scatt}}\$, \$N_{\\text{bb}}\$"
#    mod_vals+=" &  &  & 0.2 &  & &  &  & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Fscatt and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$E_{\\text{line}}\$"
#    mod_vals+=" &  &  & 0.2 &  & &  &  & & & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## Fscatt and sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\sigma_{\\text{line}}\$"
#	mod_vals+=" &  &  & 0.2 &  & &  &  & &  & & 0.01 "  ## Fscatt and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$N_{\\text{line}}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 &  &  & &  & .6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## Gamma and Tbb
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 &  &  & &  &  & 2500 &  & &  "  ## Gamma and Nbb
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & &  "  ## Fscatt and Gamma and Tbb
#    varpar="\$F_\\text{scatt}\$,  \$\\Gamma\$, \$T_\\text{bb}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & &  & 2500 &  & &  "  ## Fscatt and Gamma and Nbb
#    varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{bb}}\$ "
#	mod_vals+=" &  & & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & 2000 &  & &  "  ## Fscatt and Tbb and Nbb
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{bb}}\$, \$N_{\\text{bb}}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & &  &  &  & & 0.01 "  ## Gamma and Fscatt and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$"
#	mod_vals+=" &  & & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0  &  &  & & 0.01 "  ## Fscatt and Tbb and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{bb}}\$, \$N_{\\text{line}}\$"
#    mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & &  &  & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Fscatt and Gamma and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$E_{\\text{line}}\$"
#    mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & &  &  & & & 0.7 .005 0.1 0.1 1.0 1.0  &  "  ## Fscatt and Gamma and sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$\\sigma_{\\text{line}}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & &  &  & &  & & 0.01 "  ## Fscatt and Gamma and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$"
#    mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & 2000 &  & &  "  ## Gamma and Fscatt and Tbb and Nbb
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{bb}}\$, \$N_{\\text{bb}}\$"
#    mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Gamma and Fscatt and Tbb and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{bb}}\$, \$E_{\\text{line}}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & & & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## Gamma and Fscatt and Tbb and sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{bb}}\$, \$\\sigma_{\\text{line}}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & & & & 0.01 "  ## Gamma and Fscatt and Tbb and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$T_{\\text{bb}}\$, \$N_{\\text{line}}\$"



## For phabs*(simpler*diskpbb+gauss)
#	mod_vals+=" &  &  &  &  & & & & & &  "  ## All tied
#   varpar=" - "
#	mod_vals+=" &  &  & 0.2 &  &  & &  &  &  &  "   ## Fscatt
#	varpar="\$F_{\\text{scatt}}\$"
#	mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & &  & &  &  &  & &  "  ## Gamma
#	varpar="\$\\Gamma\$"
#	mod_vals+=" &  &  &  &  & 0.6 .002 0.1 0.1 1.0 1.0 & &  &  & & " ## TdiskP
#	varpar="\$T_{\\text{diskP}}\$"
#	mod_vals+=" &  &  &  &  &  & 0.75 & & &   & " ## p
#	varpar="\$p\$"
#	mod_vals+=" &  &  &  &  &  & & 3000 & &  & " ## NdiskP
#	varpar="\$N_{\\text{diskP}}\$"
#	mod_vals+=" &  & &  &  & &  &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Eline
#	varpar="\$E_{\\text{line}}\$"
#	mod_vals+=" &  & &  &  & &  &  &  & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## sigma
#	varpar="\\textsc{gauss} \$\\sigma_{\\text{line}}\$"
#   mod_vals+=" &  & &  &  & &  &  & &  & 0.01 "  ## Nline
#	varpar="\$N_{\\text{line}}\$"
#   mod_vals+=" &  & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & & & &  &  &  "  ## Fscatt and Gamma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$"
#	mod_vals+=" & &  & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & &  & & " ## Fscatt and TdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$"
#	mod_vals+=" &  &  & 0.2 & & & 0.75 & & & & " ## Fscatt and p
#	varpar="\$F_{\\text{scatt}}\$, \$p\$"
#	mod_vals+=" &  &  & 0.2 & &  & & 2000 &  & & " ## Fscatt and NdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$N_{\\text{diskP}}\$"
#   mod_vals+=" &  &  & 0.2 & & & &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Fscatt and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$E_{\\text{line}}\$"
#	mod_vals+="& & & 0.2 & & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & " ## Fscatt and sigma
#	varpar="\$F_{\\text{scatt}}\$, \\textsc{gauss} \$\\sigma_{\\text{line}}\$"
#	mod_vals+=" &  &  & 0.2 & & & &  &  & & 0.01 "  ## Fscatt and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$N_{\\text{line}}\$"
#	mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & &  & & " ## Fscatt and Gamma and TdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \\textsc{diskpbb} \$T_{\\text{diskP}}\$"
   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & & .75 & &  & & " ## Fscatt and Gamma and p
	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$p\$"
#   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & & & 2000 &  & & " ## Fscatt and Gamma and NdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{diskP}}\$"
#   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## Fscatt and Gamma and LineE
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$E_{\\text{line}}\$"
#   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & & & & & & 0.01 " ## Fscatt and Gamma and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$N_{\\text{line}}\$"
#   mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 &  & & " ## Fscatt and TdiskP and NdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$N_{\\text{diskP}}\$"
#   mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## Fscatt and TdiskP and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$E_{\\text{line}}\$"
#   mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & & & & 0.01" ## Fscatt and TdiskP and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$N_{\\text{line}}\$"
#   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 &  & & " ## Fscatt and TdiskP and NdiskP and Gamma
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$N_{\\text{diskP}}\$, \$\\Gamma\$"
#   mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & .75 & 2000 & & & " ## Fscatt and TdiskP and NdiskP and p
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$N_{\\text{diskP}}\$, \$p\$"
#   mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## Fscatt and TdiskP and NdiskP and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$N_{\\text{diskP}}\$, \$E_{\\text{line}}\$"
#	mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 & & & 0.01 " ## Fscatt and TdiskP and NdiskP and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$T_{\\text{diskP}}\$, \$N_{\\text{diskP}}\$, \$N_{\\text{line}}\$"
#   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0  & .75 & &  & & " ## Fscatt and Gamma and p and TdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$p\$, \$T_{\\text{diskP}}\$"
#   mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & &  & .75 & 1000 &  & & " ## Fscatt and Gamma and p and NdiskP
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$p\$, \$N_{\\text{diskP}}\$"
#	mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & &  & .75 & & 6.4 0.005 5.5 5.5 7.0 7.0  & & " ## Fscatt and Gamma and p and Eline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$p\$, \$E_{\\text{line}}\$"
#	mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & &  & .75 & & & 0.8 0.002 0.3 0.3 1.0 1.0 & " ## Fscatt and Gamma and p and sigma
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$p\$, \$\\sigma_{\\text{line}}\$"
#	mod_vals+=" & & 2.6 0.01 1.0 1.0 3.1 3.1 & 0.2 & &  & .75 & & & & 0.01" ## Fscatt and Gamma and p and Nline
#	varpar="\$F_{\\text{scatt}}\$, \$\\Gamma\$, \$p\$, \$N_{\\text{line}}\$"


##
## For phabs*(nthcomp+bbodyrad+bbodyrad+gauss)
##
#	mod_vals+=" &  &     &  &  & &  & &  &  &  & & & & "
#	mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 &  &  & &  & &  &  &  & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 &  &  &  & & & & "
#   mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 &  &  & &  & .2 &  &  &  & & & & "
#	mod_vals+=" &  &  2.6 0.01 1.0 1.0 3.1 3.1  &  &  & &  & .2 &  &  & .6 .005 0.1 0.1 2.0 2.0 & & & & "
#	mod_vals+=" &  &  2.6 0.01 1.0 1.0 3.1 3.1  &  &  & &  & &  &  & .6 .005 0.1 0.1 2.0 2.0 & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 & .6 .005 0.1 0.1 2.0 2.0 &  &  & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 &    & 3000 &  & & & & "

##
##  For P70080: PHABS*(CUTOFFPL+BBODYRAD+DISKLINE)
##
#    mod_vals+=" & &  &  & 0.18 0.005 0.1 0.1 0.3 0.3 & 0.8 0.005 0.6 0.6 1.0 1.0 & 60 0.5 40 40 100 100 &  & & & & & "
#    mod_vals+=" & &  &  & 0.18 0.005 0.1 0.1 0.3 0.3 & & 60 0.5 40 40 100 100 &  & & & & & 0.0015"
#    mod_vals+=" & &  &  & 0.18 0.005 0.1 0.1 0.3 0.3 & & 60 0.5 40 40 100 100 &  & & & & & "

##  For P70080: phabs*(nthcomp+bbodyrad+gauss)
#    mod_vals+=" & & 1.9 0.01 1.0 1.0 3.1 3.1 & & & & & & & & & & "
#    mod_vals+=" & & & & & & & .05 & & & & & "
#    mod_vals+=" & & & & & & & & .7 & & & & "
#    mod_vals+=" & & & & & & & & & 100 & & & "
#    mod_vals+=" & & & & & & & .05 & .7 & & & & "
#    mod_vals+=" & & & & & & & .05 & .7 & 100 & & & "
#    mod_vals+=" & & 1.9 0.01 1.0 1.0 3.1 3.1 & & & & & .05 & .7 & & & & "

	((i+=1))
done
((i-=1))

n_params=$(python -c "from tools import get_num_of_params; print get_num_of_params('$mod_vals', '$n_spectra')")

#echo "n params: $n_params"
export n_params

################################################################################
mcmc_file="${prefix}_${day}_${fit_specifier}_MCMC.fits"
#if [ -e "$xspec_log" ]; then rm "$xspec_log"; fi; touch "$xspec_log"
#if [ -e "${out_dir}/$mcmc_file" ]; then rm "${out_dir}/$mcmc_file"; fi;
fzpar=" - "

########################################
## Writing the rest of the xspec script
########################################

#echo "ignore 1-$i: **-3.0 25.0-**" >> $xspec_script
echo "ignore 1-$i: **-3.0 20.0-**" >> $xspec_script
#echo "notice 1-$i: 3.0-25.0" >> $xspec_script
echo "notice 1-$i: 3.0-20.0" >> $xspec_script
echo "ignore 1-$i: 11" >> $xspec_script
echo "setplot energy" >> $xspec_script
echo "systematic 0.005" >> $xspec_script
echo "xsect vern" >> $xspec_script
echo "abund wilm" >> $xspec_script

##
## GX339-4 spectral model #1:  simpler * diskbb + gauss
##
#model_string="phabs*(simpler*diskbb+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6 -0.1" >> $xspec_script ## From Reynolds and Miller 2013
#echo "newpar 2 2.6 0.005 1.0 1.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1 -0.1" >> $xspec_script
#echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
##echo "newpar 5 0.830878 -0.1" >> $xspec_script
##fzpar="\$T_{\\text{disk}}\$"
#echo "newpar 6 2505.72" >> $xspec_script
##echo "newpar 6 2505.72 -0.1" >> $xspec_script
##fzpar="\$N_{\\text{disk}}=2505.72\$"
#echo "newpar 7 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
##echo "newpar 8 0.8 .005 0.1 0.1 1.0 1.0" >> $xspec_script
#echo "newpar 8 0.97 -0.1" >> $xspec_script  ## Value from steppar on mean spectrum
#fzpar+=", \$\\sigma = 0.97\$ "
#echo "newpar 9 1.0E-02" >> $xspec_script

##
## GX339-4 spectral model #2:  simpler * diskbb + bbodyrad + gauss
##
#model_string="phabs*(simpler*diskbb+bbodyrad+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.6 0.01 1.0 1.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1" >> $xspec_script
#echo "freeze 4" >> $xspec_script
#echo "newpar 5 0.83 0.002 0.6 0.6 1.0 1.0" >> $xspec_script
###echo "newpar 5 0.955420" >> $xspec_script
##echo "newpar 5 0.830759333333 -0.1" >> $xspec_script
##fzpar="\$T_{\\text{disk}}=0.8307593\$"
#echo "newpar 6 2500 " >> $xspec_script
##echo "newpar 6 2505.72 -1" >> $xspec_script
##fzpar="\$N_{\\text{disk}} = 2505.72\$"
#echo "newpar 7 0.6 0.002 0.1 0.1 1.0 1.0" >> $xspec_script
##echo "newpar 7 0.537252 -0.1" >> $xspec_script
##fzpar="\$T_{\\text{bb}}\$"
##echo "newpar 8 600 5 50 50 1000 1000" >> $xspec_script
##echo "newpar 8 1000" >> $xspec_script
#echo "newpar 8 8857.68 -0.1" >> $xspec_script
#fzpar="\$N_{\\text{bb}}=8857.68\$"
#echo "newpar 9 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
##echo "newpar 10 0.7 .005 0.1 0.1 1.0 1.0" >> $xspec_script
#echo "newpar 10 0.82 -0.1" >> $xspec_script  ## Value from steppar on mean spectrum
##echo "newpar 10 0.56 -0.1" >> $xspec_script  ## Value from steppar on mean spectrum
#fzpar+=", \$\\sigma_\\text{line}=0.82\$"
#echo "newpar 11 0.01" >> $xspec_script

##
## GX339-4 spectral model #3: simpler * diskPbb + gauss
##
#model_string="phabs*(simpler*diskpbb+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6 -0.1" >> $xspec_script
#echo "newpar 2 2.6 0.01 1.0 1.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1 -0.1" >> $xspec_script
#echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
##echo "newpar 5 0.744958 -0.1" >> $xspec_script
##echo "newpar 5 0.838785 -0.1" >> $xspec_script
##fzpar="\$T_\\text{diskP}\$"
#echo "newpar 6 0.75" >> $xspec_script
##echo "newpar 6 0.75 -0.1" >> $xspec_script
##fzpar="\$ p = 0.75\$"
##echo "newpar 7 560.907" >> $xspec_script
#echo "newpar 7 560.907 -0.1" >> $xspec_script
#fzpar="\$N_{\\text{diskP}}=560.907\$"
#echo "newpar 8 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
##echo "newpar 9 0.7 .005 0.1 0.1 1.0 1.0" >> $xspec_script
#echo "newpar 9 0.76 -0.1" >> $xspec_script  ## Value from steppar on mean spectrum
#fzpar+=", \$\\sigma_\\text{line}=0.76\$"
#echo "newpar 10 0.01" >> $xspec_script

##
## FOR GX339-4 2010 outburst, using bbodyrad
##
#model_string="phabs*(simpler*bbodyrad+bbodyrad+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script  ## This value is from Reynolds and Miller 2013
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.6 0.01 1.0 1.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1" >> $xspec_script
#echo "freeze 4" >> $xspec_script
##echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
#echo "newpar 5 0.791511" >> $xspec_script
##echo "newpar 5 0.787971" >> $xspec_script
#echo "freeze 5" >> $xspec_script
##echo "newpar 6 3705.67" >> $xspec_script
#echo "newpar 6 1909.83" >> $xspec_script
##echo "freeze 6" >> $xspec_script
#echo "newpar 7 0.6 0.002 0.1 0.1 1.0 1.0" >> $xspec_script
##echo "newpar 8 41716.1" >> $xspec_script
#echo "newpar 8 24310.8" >> $xspec_script
#echo "freeze 8" >> $xspec_script
#echo "newpar 9 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
#echo "newpar 10 0.5 .005 0.1 0.1 0.8 0.8" >> $xspec_script
#echo "newpar 11 0.01" >> $xspec_script


##
## FOR SAXJ1808, TONY'S MODEL
##
#model_string="phabs*(cutoffpl+bbodyrad+diskline)"
#echo "mod $model_string $mod_vals" >> $xspec_script
#echo "newpar 1 0.113" >> $xspec_script  # From Wilkinson et al 2011 from Ibragimov and Poutanen 2009
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 1.96 0.005 1.1 1.1 2.1 2.1" >> $xspec_script
#echo "newpar 3 500" >> $xspec_script
##echo "freeze 3" >> $xspec_script
#echo "newpar 5 0.8 0.005 0.6 0.6 1.0 1.0" >> $xspec_script
#echo "newpar 6 60 0.5 40 40 100 100" >> $xspec_script
#echo "newpar 7 6.5 .01 6.4 6.4 6.9 6.9" >> $xspec_script
#echo "newpar 8 -3.3" >> $xspec_script
#echo "newpar 9 11.9" >> $xspec_script
#echo "newpar 10 1000" >> $xspec_script
#echo "newpar 11 60" >> $xspec_script
#echo "freeze 11" >> $xspec_script
#echo "newpar 12 0.0015" >> $xspec_script

#model_string="phabs*(nthcomp+bbodyrad+bbodyrad+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script  ## This value is from Reynolds and Miller 2013
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.6 0.01 1.0 1.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 80.0 50.0 50.0 200.0 200.0" >> $xspec_script
#echo "thaw 3" >> $xspec_script
#echo "newpar 4 0.8" >> $xspec_script
#echo "thaw 4" >> $xspec_script
#echo "newpar 5 0" >> $xspec_script
#echo "freeze 5" >> $xspec_script
#echo "newpar 6 0.0" >> $xspec_script
#echo "freeze 6" >> $xspec_script
#echo "newpar 7 0.5" >> $xspec_script
#echo "newpar 8 .8 .005 0.5 0.5 1.0 1.0" >> $xspec_script
#echo "newpar 9 3300 " >> $xspec_script
##echo "freeze 9" >> $xspec_script
#echo "newpar 10 .6 .005 0.1 0.1 2.0 2.0" >> $xspec_script
#echo "newpar 11 1.14574E+04" >> $xspec_script
#echo "freeze 11" >> $xspec_script
##echo "newpar 11 3340.2592 " >> $xspec_script
#echo "newpar 12 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
#echo "newpar 13 0.5 .005 0.1 0.1 0.8 0.8" >> $xspec_script
#echo "newpar 14 1.0E-02" >> $xspec_script
### Tieing disk-component blackbodies together
#for (( x=1; x<=i*numpar; x++ )); do
#    n=$((x%numpar))
#    y=$((x+4))
#    if (( n == 4 )) && (( n != 0 )); then
#        echo "newpar $x = $y" >> $xspec_script
#    fi
#done

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
#echo "newpar 12 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
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

#model_string="phabs*(nthcomp+bbodyrad+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.113" >> $xspec_script ## This value is from Wilkinson & Uttley 2011; Ibragimov & Poutanen 2009
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 1.9 0.01 1.0 1.0 3.1 3.1" >> $xspec_script
##echo "newpar 3 45.0 20.0 20.0 150.0 150.0" >> $xspec_script
##echo "thaw 3" >> $xspec_script
#echo "newpar 3 45.53"  >> $xspec_script
#echo "freeze 3" >> $xspec_script
#echo "newpar 4 0.8" >> $xspec_script
#echo "thaw 4" >> $xspec_script
#echo "newpar 5 0" >> $xspec_script
#echo "freeze 5" >> $xspec_script
#echo "newpar 6 0.0" >> $xspec_script
#echo "freeze 6" >> $xspec_script
#echo "newpar 7 0.05" >> $xspec_script
#echo "newpar 8 .8 .005 0.1 0.1 3.0 3.0" >> $xspec_script
#echo "newpar 9 3300 " >> $xspec_script
#echo "newpar 10 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
#echo "newpar 11 0.5 .005 0.1 0.1 0.8 0.8" >> $xspec_script
#echo "newpar 12 1.0E-04" >> $xspec_script


echo "log $xspec_log" >> $xspec_script
echo "query yes" >> $xspec_script
echo "chatter 4" >> $xspec_script
echo "fit 1000" >> $xspec_script

####################################################
## Uncomment the following to do MCMC error finding
####################################################
if (( mcmc_flag == 1 )); then
    echo "chain type gw " >> $xspec_script
    echo "chain burn 2000" >> $xspec_script
    echo "chain walkers 1000" >> $xspec_script
    echo "chain length 100000" >> $xspec_script
    echo "chain run $mcmc_file">> $xspec_script
    total_par=$((n_params*n_spectra))
    echo "total pars: $total_par"
    echo "error maximum 10000. 2.706 1-$total_par" >> $xspec_script
fi


## Getting the chisquared and degrees of freedom from the output
echo "tclout stat" >> $xspec_script
echo "scan \$xspec_tclout \"%f\" chistat" >> $xspec_script
echo "tclout dof" >> $xspec_script
echo "scan \$xspec_tclout \"%d\" degoff" >> $xspec_script
echo "echo \$chistat >> chi.txt" >> $xspec_script
echo "echo \$degoff >> dof.txt" >> $xspec_script

## Plotting the ratio of all spectra with the model
echo "cpd /xw" >> $xspec_script
echo "setplot delete all" >> $xspec_script
echo "setplot command hardcopy ${ratio_plot}.eps/cps" >> $xspec_script
echo "plot ratio" >> $xspec_script

## Plotting the first phase spectrum and the model components
echo "chatter 0" >> $xspec_script
for (( k=i; k>=2; k-- )); do
    echo "data $k none/" >> $xspec_script
done
echo "setplot delete all" >> $xspec_script
echo "setplot command @first_spectrum_fit.pco ${model_string} ${first_plot}" >> $xspec_script
echo "plot eeufspec" >> $xspec_script
echo "exit" >> $xspec_script

cd "$out_dir"
##xspec - "$xspec_script"  ## use this one for tcl
##xspec < "$xspec_script"  ## use this one for plotting

#xspec "$xspec_script" > $dump_file

# > $dump_file  ## stick this on the end to put all output into a dump file

##tail -n 10 "$xspec_log"

#open -a "TextWrangler" "$xspec_script"
#open "${ratio_plot}.eps"
#open "${first_plot}.eps"
#open "$xspec_log"

chis=$(tail -n 1 chi.txt)
dof=$(tail -n 1 dof.txt)
echo CHISQUARED = "$chis"
echo DOF = "$dof"
#
#echo " & \\textsc{$model_string} & \$$chis / $dof\$ & $varpar & " >> $tex_tab_file
echo " & \\textsc{$model_string} & \$$chis / $dof\$ & $varpar & $fzpar \\\\ " >> $tex_tab_file
## The lag chisquared, \ref, and \\ are added in fake_qpo_spectra.py
echo ""
echo "Finished running sed_fitting.sh"
cd ..

echo "python $exe_dir/multifit_plots.py $out_dir/$xspec_log" \
        --mod_string "\"${model_string}\"" \
        -w "$parfit_file"
python $exe_dir/multifit_plots.py "$out_dir/$xspec_log" \
        --mod_string "\"${model_string}\"" \
        -w "$parfit_file"

#open "$tex_tab_file"
open "$parfit_file"
export parfit_file

## In directory 'Scripts', can run log_to_textable.ipynb to put an xspec log
## file into a LaTeX table

################################################################################