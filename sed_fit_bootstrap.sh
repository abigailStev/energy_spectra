#!/bin/bash

################################################################################
## 
## Use in rxte_reduce/ccf_bootstrap.sh with ccf_bootstrap.sh and
## simulate_qpo_bootstrap.sh.
##
## Bash script for phase-resolved spectroscopy: run energyspec.py to make phase-
## resolved energy spectra, make an XSPEC SED fitting script, run
## the script, read off fit data from log file with multifit_plots.py, and make
## plots of SED parameters changing with QPO phase, fit a function to those
## parameter variations.
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
if (( $# != 6 )); then
    echo -e "\tUsage: ./sed_fitting.sh <prefix> <dt multiple> <num "\
        "seconds> <testing> <date> <num of bootstrap>\n"
    exit
fi

## ./sed_fit_bootstrap.sh GX339-BQPO 64 64 0 150526
prefix=$1
dt=$2
numsec=$3
testing=$4
day=$5
boot_num=$6

################################################################################

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. $HEADAS/headas-init.sh
fi

home_dir=$(ls -d ~)
ccf_dir="$home_dir/Dropbox/Research/cross_correlation/out_ccf/${prefix}/bootstrapped"
exe_dir="$home_dir/Dropbox/Research/energy_spectra"
out_dir="$exe_dir/out_es/${prefix}/bootstrapped"
dump_file=dump.txt # Name of dumping file for intermediary steps
data_dir="$home_dir/Reduced_data/$prefix"
bkgd_spec="$data_dir/evt_bkgd_rebinned.pha"
rsp_matrix="${prefix}_PCU2.rsp"

spec_type=0  # 0 for mean+ccf, 1 for ccf, 2 for mean

fit_specifier+="1BB-FS-G-Tin"
#fit_specifier+="1BB-FS-G"
#fit_specifier+="pBB-FS-Tin-G"
#fit_specifier+="2BB-FS-G-kT"
fit_specifier+="-fzs-fzNbb"
#fit_specifier+="-fzs"
#fit_specifier+="-fzNbb"
#fit_specifier+="-wMCMC"
export fit_specifier

tex_tab_file="$home_dir/Dropbox/Research/CCF_paper1/ensp_models.txt"
parfit_file="$out_dir/${prefix}_${day}_${fit_specifier}_sines.txt"
tab_ext="dat"
plot_ext="eps"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ -e "$parfit_file" ]; then rm "$parfit_file"; fi; touch "$parfit_file"
if [ ! -e "$tex_tab_file" ]; then touch "$tex_tab_file"; fi

if [ -e "$data_dir/PCU2.rsp" ]; then
    cp "$data_dir/PCU2.rsp" "${out_dir}/$rsp_matrix"
else
    echo "ERROR: Response matrix doesn't exist in the reduced data directory."
    exit
fi

ccf_file="$ccf_dir/${prefix}_${day}_t${dt}_${numsec}sec_adj_b-1.fits"
if (( $testing==1 )); then
    ccf_file="$ccf_dir/test_${prefix}_${day}_t${dt}_${numsec}sec_adj_b-1.fits"
fi
obs_time=$(python -c "from tools import get_key_val; print get_key_val('$ccf_file', 0, 'EXPOSURE')")
#echo "$obs_time"

############################################
## Looping through the bootstrap iterations
############################################
for (( b=1; b<=boot_num; b++ )); do
	if (( b % 1 == 0 )); then echo -e "\t $b"; fi

	boot_fit_specifier="${fit_specifier}_b-${b}"
    out_name="${prefix}_${day}_t${dt}_${numsec}sec_adj_b-${b}"
    xspec_script="$out_dir/${prefix}_${day}_${boot_fit_specifier}_xspec.xcm"
    xspec_log="${prefix}_${day}_${boot_fit_specifier}_xspec.log"
    ratio_plot="${prefix}_${day}_${boot_fit_specifier}_ratio"
    first_plot="${prefix}_${day}_${boot_fit_specifier}_firstspec"
    spectrum_plot="${prefix}_${day}_${boot_fit_specifier}_allspectra"
    ccf_file="$ccf_dir/${out_name}.fits"
    if (( $testing==1 )); then
        ccf_file="$ccf_dir/test_${out_name}.fits"
    fi
    if [ -e "$xspec_script" ]; then rm "$xspec_script"; fi; touch "$xspec_script"

    echo "lmod simpler /Users/abigailstevens/Dropbox/Research/xspecmods/" >> $xspec_script

    i=1
    mod_vals=""

    cd "$out_dir"

    ###############################################################
    ## Generating a spectral energy distribution at each phase bin
    ###############################################################

    n_spectra=$(( 8205-8182+1 ))
    #echo "$n_spectra"
    export n_spectra

    for (( timebin=8182; timebin<=8205; timebin++ )); do
        if (( timebin>=8192 )); then
            tbin=$((timebin-8192))
        else
            tbin=$((timebin))
        fi

        out_end="${out_name}_ccfwmean_${tbin}bin"
        if [ -e "${ccf_file}" ]; then
            python "$exe_dir"/energyspec.py "${ccf_file}" \
                    "${out_end}.${tab_ext}" -b "$tbin" -s "$spec_type"
    #        echo "${ccf_file}"
    #        echo "${out_end}.${tab_ext}"
    #        open -a "TextWrangler" "${out_end}.${tab_ext}"
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
        ## Don't need to include background file since I've subtracted
        ## background count rate per energy channel from the mean count rate
        ## per CI

        else
            echo -e "\tERROR: ASCII2PHA was not run. Spectrum and/or response "\
                    "matrix do not exist."
        fi

        if [ ! -e "${out_end}.pha" ]; then
            echo -e "\tERROR: ASCII2PHA failed to create ${out_end}.pha."
            echo -e "\tExiting script."
            exit
        fi

        ## Writing file names to xspec script
        echo "data $i:$i $out_end.pha" >> $xspec_script

        #######################################################################
        ## Go through and uncomment the two lines for the untied parameter
        ## combination you want. First line is what's fed into XSPEC
        ## ("mod_vals"), second line is what's written in the table ("varpar").
        #######################################################################

        ##
        ## For phabs*(simpler*diskbb+gauss)
        ##
    #	mod_vals+="& & & & & & & & & "  ##  all tied
    #	varpar=" - "
    #	mod_vals+="& & & 0.2 & & & & & & "  ## FracSctr
    #	varpar="\\texttt{simpler} FracSctr"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & & & & & & & " ## Gamma
    #	varpar="\\texttt{simpler} Gamma"
    #	mod_vals+="& & & & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & " ## Tin
    #	varpar="\\texttt{diskbb} T\$_{\\text{in}}"
    #	mod_vals+="& & & & & & 3000 & & & " ##  norm(BB)
    #	varpar="\\texttt{diskbb} norm"
    #	mod_vals+="& & & & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## lineE
    #	varpar="\\texttt{gauss} LineE"
    #	mod_vals+="& & & & & & & & & 0.01"  ## norm(E)
    #	varpar="\\texttt{gauss} norm"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & & & " ## FracSctr and Gamma
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma"
    #	mod_vals+="& & & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & " ## FracSctr and Tin
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{diskbb} T\$_{\\text{in}}"
    #	mod_vals+="& & & 0.2 & & & 3000 & & & " ## FracSctr and norm(BB)
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{diskbb} norm"
    #	mod_vals+="& & & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## FracSctr and lineE
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} LineE"
    #	mod_vals+="& & & 0.2 & & & & & & 0.01"  ## FracSctr and norm(E)
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} norm"
        mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & "  ## FracSctr, Gamma, Tin
        varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} T\$_{\\text{in}}\$"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & 3000 & & & "  ## FracSctr, Gamma, norm(BB)
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{diskbb} norm"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & " ## FracSctr, Gamma, lineE
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} LineE"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & & & 0.01"  ## FracSctr, Gamma, norm(E)
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & & & & 0.01 "  ## Gamma, FracSctr, norm(E), Tin
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{diskbb} T\$_{\\text{in}}\$"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & 0.8 0.002 0.3 0.3 1.0 1.0 & 3000 & & & "  ## Gamma, FracSctr, norm(E), and norm(BB)
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{diskbb} norm"
    #	mod_vals+="& & 2.6 0.005 2.0 2.0 3.1 3.1 & 0.2 & & & & 6.4 0.005 5.5 5.5 7.0 7.0 & & 0.01 " ## FracSctr, Gamma, norm(E), lineE
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm, \\texttt{gauss} LineE"

        ##
        ## For phabs*(simpler*diskbb+bbodyrad+gauss) or phabs*(simpler*bbodyrad+bbodyrad+gauss)
        ##
    #	mod_vals+=" &  &  &     &  & &  &  &  &  & &  "  ## All tied
    #	varpar=" - "
    #	mod_vals+=" &  &  & 0.2 &  & &  &  &  &  & &   "   ## FracSctr
    #	varpar="\\texttt{simpler} FracSctr"
    #	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & &  & &  &  &  &  & &  "  ## Gamma
    #	varpar="\\texttt{simpler} Gamma"
    #	mod_vals+=" &  &  &  &  & &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## bb kT
    #	varpar="\\texttt{bbodyrad} kT"
    #	mod_vals+=" &  & &  &  & &  &  & 400 &  & &  "  ## bb norm
    #	varpar="\\texttt{bbodyrad} norm"
    #	mod_vals+=" &  & &  &  & &  &  &  & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## line E
    #	varpar="\\texttt{gauss} LineE"
    #    mod_vals+=" &  & &  &  & &  &  &  & &  & 0.01 "  ## E norm
    #	varpar="\\texttt{gauss} norm"
    #    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & &  &  &  & &  "  ## FracSctr and Gamma
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma"
    #	mod_vals+=" &  &  & 0.2 &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  &  &  & & " ## FracSctr and diskbb Tin
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{diskbb} T\$_{\\text{in}}\$"
    #	mod_vals+=" &  &  & 0.2 &  & &  & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## FracSctr and bb kT
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT"
    #	mod_vals+=" &  &  & 0.2 &  & &  &  & 2500 &  & &  "  ## FracSctr and bb norm
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} norm"
    #    mod_vals+=" &  &  & 0.2 &  & &  &  & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## FracSctr and line E
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} LineE"
    #	mod_vals+=" &  &  & 0.2 &  & &  &  & &  & & 0.01 "  ## FracSctr and E norm
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{gauss} norm"
    #	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 &  &  & &  & .6 .002 0.1 0.1 1.0 1.0 &  &  & & " ## Gamma and bb kT
    #	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 &  &  & &  &  & 2500 &  & &  "  ## Gamma and bb norm
    #	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 &  &  & &  "  ## Gamma and FracSctr and bb kT
    #    varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT, \\texttt{simpler} Gamma"
    #	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & &  & 2000 &  & &  "  ## Gamma and FracSctr and bb norm
    #    varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{bbodyrad} norm "
    #	mod_vals+=" &  & & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & 2000 &  & &  "  ## FracSctr and bb kT and bb norm
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT, \\texttt{bbodyrad} norm"
    #	mod_vals+=" &  & 2.8 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & &  &  &  & & 0.01 "  ## Gamma and FracSctr and E norm
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm"
    #	mod_vals+=" &  & & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0  &  &  & & 0.01 "  ## FracSctr and bb kT and E norm
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{bbodyrad} kT, \\texttt{gauss} norm"
    #    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & &  &  & & 6.4 0.005 5.5 5.5 7.0 7.0 & &  "  ## FracSctr and Gamma and line E
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} LineE"
    #    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & &  &  & & & 0.8 0.005 0.4 0.4 1.0 1.0 &  "  ## FracSctr and Gamma and sigma
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} Sigma"
    #	mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & &  &  & &  & & 0.01 "  ## FracSctr and Gamma and E norm
    #	varpar="\\texttt{simpler} FracSctr, \\texttt{simpler} Gamma, \\texttt{gauss} norm"
    #    mod_vals+=" &  & 2.6 0.01 2.0 2.0 3.1 3.1 & 0.2 &  & & & 0.6 .002 0.1 0.1 1.0 1.0 & 20000 &  & &  "  ## Gamma and FracSctr and bb kT and bb norm


        ##
        ## For phabs*(simpler*diskpbb+gauss)
        ##
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

        ((i+=1))
    done
    ((i-=1))

    n_params=$(python -c "from tools import get_num_of_params; print get_num_of_params('$mod_vals', '$n_spectra')")
    echo "n params: $n_params"
    export n_params

    ############################################################################
    mcmc_file="${prefix}_${day}_${boot_fit_specifier}_MCMC.fits"
    if [ -e "$xspec_log" ]; then rm "$xspec_log"; fi; touch "$xspec_log"
    #if [ -e "${out_dir}/$mcmc_file" ]; then rm "${out_dir}/$mcmc_file"; fi;
    fzpar=" "

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
    ## Simple test for GX339-4 2010 outburst
    ##
    model_string="phabs*(simpler*diskbb+gauss)"
    echo "mod ${model_string} $mod_vals" >> $xspec_script
    echo "newpar 1 0.6" >> $xspec_script ## From Reynolds and Miller 2013
    echo "freeze 1" >> $xspec_script
    echo "newpar 2 2.6 0.005 2.0 2.0 3.1 3.1" >> $xspec_script
    echo "newpar 3 0.2" >> $xspec_script
    echo "newpar 4 1" >> $xspec_script
    echo "freeze 4" >> $xspec_script
    echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
    #echo "newpar 5 0.830878" >> $xspec_script
    #echo "freeze 5" >> $xspec_script
    #fzpar="\\texttt{diskbb} T\$_{\\text{in}}\$"
    echo "newpar 6 2505.72" >> $xspec_script
    echo "freeze 6" >> $xspec_script
    fzpar="\\texttt{diskbb} norm=2505.72"
    echo "newpar 7 6.4 0.005 5.5 5.5 7.0 7.0" >> $xspec_script
    #echo "newpar 8 0.5 .005 0.1 0.1 0.8 0.8" >> $xspec_script
    echo "newpar 8 0.97" >> $xspec_script  ## Value from steppar on mean spectrum
    echo "freeze 8" >> $xspec_script
    fzpar+="\\texttt{gauss} sigma=0.97 "
    echo "newpar 9 1.0E-02" >> $xspec_script

    ##
    ## FOR GX339-4 2010 outburst, using diskBB
    ##
    #model_string="phabs*(simpler*diskbb+bbodyrad+gauss)"
    #echo "mod ${model_string} $mod_vals" >> $xspec_script
    #echo "newpar 1 0.6" >> $xspec_script
    #echo "freeze 1" >> $xspec_script
    #echo "newpar 2 2.6 0.01 2.0 2.0 3.1 3.1" >> $xspec_script
    #echo "newpar 3 0.2" >> $xspec_script
    #echo "newpar 4 1" >> $xspec_script
    #echo "freeze 4" >> $xspec_script
    ##echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
    ##echo "newpar 5 0.744958" >> $xspec_script
    #echo "newpar 5 0.828058" >> $xspec_script
    ##echo "freeze 5" >> $xspec_script
    ##echo "newpar 6 3214.06" >> $xspec_script
    #echo "newpar 6 2519.16" >> $xspec_script
    ##echo "freeze 6" >> $xspec_script
    #echo "newpar 7 0.6 0.002 0.1 0.1 1.0 1.0" >> $xspec_script
    #echo "newpar 8 2500" >> $xspec_script
    ##echo "newpar 8 473.804" >> $xspec_script
    ##echo "newpar 8 671.474" >> $xspec_script
    ##echo "freeze 8" >> $xspec_script
    #echo "newpar 9 6.5 0.002 6.4 6.4 6.8 6.8" >> $xspec_script
    #echo "newpar 10 0.7 .005 0.1 0.1 1.0 1.0" >> $xspec_script
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
    ##echo "newpar 5 0.8 0.002 0.5 0.5 1.0 1.0" >> $xspec_script
    ##echo "newpar 5 0.744958" >> $xspec_script
    #echo "newpar 5 0.838785" >> $xspec_script
    #echo "freeze 5" >> $xspec_script
    #echo "newpar 6 0.75" >> $xspec_script
    ##echo "newpar 7 3214.06" >> $xspec_script
    ##echo "freeze 7" >> $xspec_script
    #echo "newpar 8 6.5 0.002 6.4 6.4 6.8 6.8" >> $xspec_script
    #echo "newpar 9 0.7 .005 0.1 0.1 1.0 1.0" >> $xspec_script
    #echo "newpar 10 0.01" >> $xspec_script

    ##
    ## The rest of the normal commands, and telling it to fit!
    ##
    echo "log $xspec_log" >> $xspec_script
    echo "query yes" >> $xspec_script
    echo "chatter 4" >> $xspec_script
    echo "fit 1000" >> $xspec_script

    ##
    ## Uncomment the following to do MCMC error finding
    ##
    #echo "chain type gw " >> $xspec_script
    #echo "chain burn 2000" >> $xspec_script
    #echo "chain walkers 1000" >> $xspec_script
    #echo "chain length 100000" >> $xspec_script
    #echo "chain run $mcmc_file">> $xspec_script
    #echo "error maximum 10000. 2.706 1-216" >> $xspec_script
    ##echo "error maximum 10000. 2.706 1-240" >> $xspec_script
    ##echo "error maximum 10000. 2.706 1-264" >> $xspec_script
    ##echo "error maximum 10000. 2.706 1-288" >> $xspec_script

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
    xspec "$xspec_script" > $dump_file  ## stick this on the end to put all output into a dump file

    ##tail -n 10 "$xspec_log"
    #open -a "TextWrangler" "$xspec_script"
    #open "${ratio_plot}.eps"
    #open "${first_plot}.eps"
    #open "$xspec_log"

    chis=$(tail -n 1 chi.txt)
    dof=$(tail -n 1 dof.txt)
    #echo CHISQUARED = "$chis"
    #echo DOF = "$dof"

    echo " & \\texttt{$model_string} & \$$chis / $dof\$ & $varpar & $fzpar & " >> $tex_tab_file
    ## The lag chisquared, \ref, and \\ are added in fake_qpo_spectra.py
    #echo ""
    #echo "Finished running sed_fitting.sh"
    cd ..

    #############################
    ## Running multifit_plots.py
    #############################

    python $exe_dir/multifit_plots.py "$out_dir/$xspec_log" \
            --mod_string "\"${model_string}\"" \
            -w "$parfit_file" \
            --quiet

done

#open "$tex_tab_file"
open "$parfit_file"
export parfit_file

## In directory 'Scripts', can run log_to_textable.ipynb to put an xspec log
## file into a LaTeX table

################################################################################
