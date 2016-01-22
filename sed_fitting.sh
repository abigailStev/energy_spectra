#!/bin/bash

################################################################################
## 
## Bash script for phase-resolved spectroscopy: run energyspec.py to make phase-
## resolved energy spectra, make an XSPEC SED fitting script, run
## the script, read off fit data from log file with multispec_plots.py, and make
## plots of SED parameters changing with QPO phase, fit a function to those
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
## Written by Abigail Stevens <A.L.Stevens at uva.nl> 2015-2016
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

spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean
mcmc_flag=0  # 0 for no, 1 for yes

fit_specifier="1BB-FS-G-NE"
#fit_specifier="pBB-FS-Tin-G"
#fit_specifier="2BB-FS-G-Nbb"
#fit_specifier="2BB-FS-G-kT"
#fit_specifier+="-fzs-fzNbb"
#fit_specifier+="-fzs"
#fit_specifier+="-fzNbb"
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
tex_tab_file="$home_dir/Dropbox/Research/CCF_paper1/spec_models_1BB.txt"
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
if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"
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
#n_spectra=47
#echo "$n_spectra"
export n_spectra

for (( timebin=8182; timebin<=8205; timebin++ )); do
#for (( timebin=8182; timebin<=8228; timebin++ )); do

    if (( timebin>=8192 )); then
        tbin=$((timebin-8192))
    else
        tbin=$((timebin))
    fi

	out_end="${out_name}_ccfwmean_${tbin}bin"
    if [ ! -e "${out_end}.pha" ]; then
        if [ -e "${ccf_file}" ]; then
            python "$exe_dir"/energyspec.py "${ccf_file}" \
                    "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
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
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & & & & & & & " ## Gamma
#	varpar="\\texttt{simpler} Gamma"
#	mod_vals+="& & & 0.2 & & & & & & "  ## FracSctr
#	varpar="\\texttt{simpler} FracSctr"
#	mod_vals+="& & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & " ## Tin
#	varpar="\\texttt{diskbb} T\$_{\\text{in}}\$"
#	mod_vals+="& & & & & & 3000 & & & " ##  norm(BB)
#	varpar="\\texttt{diskbb} norm"
#	mod_vals+="& & & & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## lineE
#	varpar="\\texttt{gauss} LineE"
#	mod_vals+="& & & & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & " ## sigma
#	varpar="\\texttt{gauss} sigma"
#	mod_vals+="& & & & & & & & & 0.01"  ## norm(E)
#	varpar="\\texttt{gauss} norm"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & & & " ## FracSctr and Gamma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma"
#	mod_vals+="& & & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & " ## FracSctr and Tin
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskbb} T\$_{\\text{in}}\$"
#	mod_vals+="& & & 0.2 & & & 3000 & & & " ## FracSctr and norm(BB)
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskbb} norm"
#	mod_vals+="& & & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## FracSctr and lineE
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} LineE"
#	mod_vals+="& & & 0.2 & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & " ## FracSctr and sigma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} sigma"
#	mod_vals+="& & & 0.2 & & & & & & 0.01"  ## FracSctr and norm(E)
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} norm"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & "  ## FracSctr, Gamma, Tin
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} T\$_{\\text{in}}\$"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & 3000 & & & "  ## FracSctr, Gamma, norm(BB)
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} norm"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## FracSctr, Gamma, lineE
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} LineE"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & "  ## FracSctr, Gamma, sigma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} sigma"
 	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & & & 0.01"  ## FracSctr, Gamma, norm(E)
	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm"
#    mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & 3000 & & & "  ## FracSctr, Gamma, Tin, norm(BB)
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} T\$_{\\text{in}}\$, \\texttt{diskbb} norm"
#    mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & 6.4 0.005 5.5 5.5 7.0 7.0 & & "  ## FracSctr, Gamma, Tin, lineE
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} T\$_{\\text{in}}\$, \\texttt{gauss} LineE"
#    mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & 0.8 0.002 0.3 0.3 1.0 1.0 & "  ## FracSctr, Gamma, Tin, lineE
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} T\$_{\\text{in}}\$, \\texttt{gauss} sigma"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & 0.01 "  ## Gamma, FracSctr, norm(E), Tin
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{diskbb} T\$_{\\text{in}}\$"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & 3000 & & & 0.1 "  ## Gamma, FracSctr, norm(E), and norm(BB)
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{diskbb} norm"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & 0.01 " ## FracSctr, Gamma, norm(E), lineE
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{gauss} LineE"
#	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & 0.01 " ## FracSctr, Gamma, norm(E), lineE
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{gauss} sigma"

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
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & &  & &  &  &  &  & &  "  ## Gamma
#	varpar="\\texttt{simpler} Gamma"
#	mod_vals+=" &  &  & 0.2 &  & &  &  &  &  & &   "   ## FracSctr
#	varpar="\\texttt{simpler} FracSctr"
#	mod_vals+=" &  &  &  &  & &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## bb kT
#	varpar="\\texttt{bbodyrad} kT"
#	mod_vals+=" &  & &  &  & &  &  & 400 &  & &  "  ## bb norm
#	varpar="\\texttt{bbodyrad} norm"
#	mod_vals+=" &  & &  &  & &  &  &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## line E
#	varpar="\\texttt{gauss} LineE"
#	mod_vals+=" &  & &  &  & &  &  &  &  & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## sigma
#	varpar="\\texttt{gauss} sigma"
#    mod_vals+=" &  & &  &  & &  &  &  & &  & 0.01 "  ## E norm
#	varpar="\\texttt{gauss} norm"
#    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & &  &  &  & &  "  ## FracSctr and Gamma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma"
#	mod_vals+=" &  &  & 0.2 &  & &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## FracSctr and bb kT
#	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT"
#	mod_vals+=" &  &  & 0.2 &  & &  &  & 2500 &  & &  "  ## FracSctr and bb norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} norm"
#    mod_vals+=" &  &  & 0.2 &  & &  &  & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## FracSctr and line E
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} LineE"
#    mod_vals+=" &  &  & 0.2 &  & &  &  & & & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## FracSctr and sigma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} sigma"
#	mod_vals+=" &  &  & 0.2 &  & &  &  & &  & & 0.01 "  ## FracSctr and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} norm"
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 &  &  & &  & .6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## Gamma and bb kT
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 &  &  & &  &  & 2500 &  & &  "  ## Gamma and bb norm
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & &  "  ## FracSctr and Gamma and bb kT
#    varpar="\\texttt{simpler} FracSctr,  \\texttt{simpler} Gamma, \\texttt{bbodyrad} kT"
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & &  & 2500 &  & &  "  ## FracSctr and Gamma and bb norm
#    varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{bbodyrad} norm "
#	mod_vals+=" &  & & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & 2000 &  & &  "  ## FracSctr and bb kT and bb norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT, \\texttt{bbodyrad} norm"
#	mod_vals+=" &  & 2.8 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & &  &  &  & & 0.01 "  ## Gamma and FracSctr and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm"
#	mod_vals+=" &  & & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0  &  &  & & 0.01 "  ## FracSctr and bb kT and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT, \\texttt{gauss} norm"
#    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & &  &  & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## FracSctr and Gamma and line E
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} LineE"
#    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & &  &  & & & 0.7 .005 0.1 0.1 1.0 1.0  &  "  ## FracSctr and Gamma and sigma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} sigma"
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & &  &  & &  & & 0.01 "  ## FracSctr and Gamma and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm"
#    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & 2000 &  & &  "  ## Gamma and FracSctr and bb kT and bb norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{bbodyrad} kT, \\texttt{bbodyrad} norm"
#    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## Gamma and FracSctr and bb kT and line E
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{bbodyrad} kT, \\texttt{gauss} LineE"
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & & & 0.7 .005 0.1 0.1 1.0 1.0 &  "  ## Gamma and FracSctr and bb kT and sigma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{bbodyrad} kT, \\texttt{gauss} sigma"
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & & & & 0.01 "  ## Gamma and FracSctr and bb kT and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{bbodyrad} kT, \\texttt{gauss} norm"



## For phabs*(simpler*diskpbb+gauss)
#	mod_vals+=" &  &  &  &  & & & & & &  "  ## All tied
#    varpar=" - "
#	mod_vals+=" &  &  & 0.2 &  &  & &  &  &  &  "   ## FracSctr
#	varpar="\\texttt{simpler} FracSctr"
#	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & &  & &  &  &  & &  "  ## Gamma
#	varpar="\\texttt{simpler} Gamma"
#	mod_vals+=" &  &  &  &  & 0.6 .002 0.1 0.1 1.0 1.0 & &  &  & & " ## diskpbb Tin
#	varpar="\\texttt{diskpbb} T\$_{\\text{in}}\$"
#	mod_vals+=" &  &  &  &  &  & 0.75 & & &   & " ## diskpbb p
#	varpar="\\texttt{diskpbb} p"
#	mod_vals+=" &  &  &  &  &  & & 3000 & &  & " ## diskpbb norm
#	varpar="\\texttt{diskpbb} norm"
#	mod_vals+=" &  & &  &  & &  &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## line E
#	varpar="\\texttt{gauss} LineE"
#    mod_vals+=" &  & &  &  & &  &  & &  & 0.01 "  ## E norm
#	varpar="\\texttt{gauss} norm"
#    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 & & & & &  &  &  "  ## FracSctr and Gamma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma"
#	mod_vals+=" & &  & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & &  & & " ## FracSctr and diskpbb Tin
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$"
#	mod_vals+=" &  &  & 0.2 & & & 0.75 & & & & " ## FracSctr and diskpbb p
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} p"
#	mod_vals+=" &  &  & 0.2 & &  & & 2000 &  & & " ## FracSctr and diskpbb norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} norm"
#    mod_vals+=" &  &  & 0.2 & & & &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## FracSctr and line E
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} LineE"
#	mod_vals+=" &  &  & 0.2 & & & &  &  & & 0.01 "  ## FracSctr and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} norm"
#	mod_vals+=" & & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & &  & & " ## FracSctr and diskpbb Tin and Gamma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{simpler} Gamma"
#    mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & .75 & &  & & " ## FracSctr and diskpbb Tin and diskpbb p
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{diskpbb} p"
#    mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 &  & & " ## FracSctr and diskpbb Tin and diskpbb norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{diskpbb} norm"
#    mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## FracSctr and diskpbb Tin and line E
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{gauss} LineE"
#    mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & & & & 0.01" ## FracSctr and diskpbb Tin and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{gauss} norm"
#    mod_vals+=" & & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 &  & & " ## FracSctr and diskpbb Tin and diskpbb norm and Gamma
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{diskpbb} norm, \\texttt{simpler} Gamma"
#    mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & .75 & 2000 & & & " ## FracSctr and diskpbb Tin and diskpbb norm and diskpbb p
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{diskpbb} norm, \\texttt{diskpbb} p"
#    mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## FracSctr and diskpbb Tin and diskpbb norm and line E
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{diskpbb} norm, \\texttt{gauss} LineE"
#	mod_vals+=" & & & 0.2 & & 0.6 .002 0.1 0.1 1.0 1.0 & & 2000 & & & 0.01 " ## FracSctr and diskpbb Tin and diskpbb norm and E norm
#	varpar="\\texttt{simpler} FracSctr, \\texttt{diskpbb} T\$_{\\text{in}}\$, \\texttt{diskpbb} norm, \\texttt{gauss} norm"

##
## For phabs*(nthcomp+bbodyrad+bbodyrad+gauss)
##
#	mod_vals+=" &  &     &  &  & &  & &  &  &  & & & & "
#	mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 &  &  & &  & &  &  &  & & & & "
#	mod_vals+=" &  &     &  &  & &  & .2 &  &  &  & & & & "
#    mod_vals+=" &  & 2.8 0.01 2.6 2.6 3.1 3.1 &  &  & &  & .2 &  &  &  & & & & "
#	mod_vals+=" &  &  2.8 0.01 2.0 2.0 3.1 3.1  &  &  & &  & .2 &  &  & .6 .005 0.1 0.1 2.0 2.0 & & & & "
#	mod_vals+=" &  &  2.8 0.01 2.6 2.6 3.1 3.1  &  &  & &  & &  &  & .6 .005 0.1 0.1 2.0 2.0 & & & & "
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
#     mod_vals+=" & & & & & & & .05 & .7 & & & & "
#     mod_vals+=" & & & & & & & .05 & .7 & 100 & & & "
#    mod_vals+=" & & 1.9 0.01 1.0 1.0 3.1 3.1 & & & & & .05 & .7 & & & & "

	((i+=1))
done
((i-=1))

n_params=$(python -c "from tools import get_num_of_params; print get_num_of_params('$mod_vals', '$n_spectra')")

#echo "n params: $n_params"
export n_params

################################################################################
mcmc_file="${prefix}_${day}_${fit_specifier}_MCMC.fits"
if [ -e "$xspec_log" ]; then rm "$xspec_log"; fi; touch "$xspec_log"
if [ -e "${out_dir}/$mcmc_file" ]; then rm "${out_dir}/$mcmc_file"; fi;
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
model_string="phabs*(simpler*diskbb+gauss)"
echo "mod ${model_string} $mod_vals" >> $xspec_script
echo "newpar 1 0.6 -0.1" >> $xspec_script ## From Reynolds and Miller 2013
echo "newpar 2 2.6 0.005 2.0 2.0 3.1 3.1" >> $xspec_script
echo "newpar 3 0.2" >> $xspec_script
echo "newpar 4 1" >> $xspec_script
echo "freeze 4" >> $xspec_script
echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
#echo "newpar 5 0.830878" >> $xspec_script
#echo "freeze 5" >> $xspec_script
#fzpar="\\texttt{diskbb} T\$_{\\text{in}}\$"
echo "newpar 6 2505.72" >> $xspec_script
#echo "freeze 6" >> $xspec_script
#fzpar="\\texttt{diskbb} norm\$\\,=2505.72\$"
echo "newpar 7 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
echo "newpar 8 0.8 .005 0.1 0.1 1.0 1.0" >> $xspec_script
#echo "newpar 8 0.97" >> $xspec_script  ## Value from steppar on mean spectrum
#echo "freeze 8" >> $xspec_script
#fzpar+=", \\texttt{gauss} \$\\sigma = 0.97\$ "
echo "newpar 9 1.0E-02" >> $xspec_script

##
## GX339-4 spectral model #2:  simpler * diskbb + bbodyrad + gauss
##
#model_string="phabs*(simpler*diskbb+bbodyrad+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.6 0.01 2.0 2.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1" >> $xspec_script
#echo "freeze 4" >> $xspec_script
#echo "newpar 5 0.83 0.002 0.6 0.6 1.0 1.0" >> $xspec_script
###echo "newpar 5 0.955420" >> $xspec_script
##echo "newpar 5 0.830759333333" >> $xspec_script
##echo "freeze 5" >> $xspec_script
##fzpar="\\texttt{diskbb} T\$_{\\text{in}}\$=0.830759333333"
#echo "newpar 6 2500" >> $xspec_script
##echo "freeze 6" >> $xspec_script
##fzpar="\\texttt{diskbb} norm = 2505.72"
#echo "newpar 7 0.6 0.002 0.1 0.1 1.0 1.0" >> $xspec_script
##echo "newpar 7 0.537252" >> $xspec_script
##echo "freeze 7" >> $xspec_script
##fzpar="\\texttt{bbodyrad} kT"
##echo "newpar 8 600 5 50 50 1000 1000" >> $xspec_script
##echo "newpar 8 1000" >> $xspec_script
#echo "newpar 8 8857.68" >> $xspec_script
#echo "freeze 8" >> $xspec_script
#fzpar="\\texttt{bbodyrad} norm\$\\,=8857.68\$"
#echo "newpar 9 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
##echo "newpar 10 0.7 .005 0.1 0.1 1.0 1.0" >> $xspec_script  ## Value from steppar on mean spectrum
#echo "newpar 10 0.82" >> $xspec_script  ## Value from steppar on mean spectrum
##echo "newpar 10 0.56" >> $xspec_script  ## Value from steppar on mean spectrum
#echo "freeze 10" >> $xspec_script
#fzpar+=", \\texttt{gauss} \$\\sigma=0.82\$"
#echo "newpar 11 0.01" >> $xspec_script

##
## FOR GX339-4 2010 outburst, using diskpbb
##
#model_string="phabs*(simpler*diskpbb+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.6 0.01 2.0 2.0 3.1 3.1" >> $xspec_script
#echo "newpar 3 0.2" >> $xspec_script
#echo "newpar 4 1" >> $xspec_script
#echo "freeze 4" >> $xspec_script
#echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
##echo "newpar 5 0.744958" >> $xspec_script
##echo "newpar 5 0.838785" >> $xspec_script
##echo "freeze 5" >> $xspec_script
##fzpar="\\texttt{diskpbb} T\$_{\\text{in}}\$"
#echo "newpar 6 0.525230" >> $xspec_script
##echo "freeze 6" >> $xspec_script
##fzpar="\\texttt{diskpbb} p"
#echo "newpar 7 685.176" >> $xspec_script
##echo "freeze 7" >> $xspec_script
##fzpar="\\texttt{diskpbb} norm"
#echo "newpar 8 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
#echo "newpar 9 0.76" >> $xspec_script  ## Value from steppar on mean spectrum
#echo "freeze 9" >> $xspec_script
#fzpar+="\\texttt{gauss} sigma"
#echo "newpar 10 0.01" >> $xspec_script

##
## FOR GX339-4 2010 outburst, using bbodyrad
##
#model_string="phabs*(simpler*bbodyrad+bbodyrad+gauss)"
#echo "mod ${model_string} $mod_vals" >> $xspec_script
#echo "newpar 1 0.6" >> $xspec_script  ## This value is from Reynolds and Miller 2013
#echo "freeze 1" >> $xspec_script
#echo "newpar 2 2.6 0.01 2.0 2.0 3.1 3.1" >> $xspec_script
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
#echo "newpar 2 2.8 0.01 2.0 2.0 3.1 3.1" >> $xspec_script
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
    total_par=$((n_params*n_spectra))
    echo "total pars: $total_par"
    echo "chain type gw " >> $xspec_script
    echo "chain burn 2000" >> $xspec_script
    echo "chain walkers 1000" >> $xspec_script
    echo "chain length 100000" >> $xspec_script
    echo "chain run $mcmc_file">> $xspec_script
    total_par=$((n_params*n_spectra))
    echo "total pars: $total_par"
    #echo "error maximum 10000. 2.706 1-$total_par" >> $xspec_script
    echo "error maximum 10000. 1 1-423" >> $xspec_script
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
xspec "$xspec_script" > $dump_file    ## use this one if you expect interactivity with fit
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
#echo " & \\texttt{$model_string} & \$$chis / $dof\$ & $varpar & " >> $tex_tab_file
echo " & \\texttt{$model_string} & \$$chis / $dof\$ & $varpar & $fzpar \\\\ " >> $tex_tab_file
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
#open "$parfit_file"
export parfit_file

## In directory 'Scripts', can run log_to_textable.ipynb to put an xspec log
## file into a LaTeX table

################################################################################