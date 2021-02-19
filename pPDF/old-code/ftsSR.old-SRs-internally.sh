#!/bin/bash 

if [ $# -ne 8 ]; then
    echo "Usage: $0 <summed ratio file> <tmin fit> <tstep fit> <tmax fit> <cfgs> <complexity (Re-1) or (Im-2)> <snk row> <src row>"
    exit 1
fi

rm thisCorr.list
rm -rf JACK/
mkdir -p JACK/res/
mkdir -p JACK/cov/

sRFile=$1 # File containing summed ratio for each tsep
tminfit=$2
tstepfit=$3
tmaxfit=$4
cfgs=$5
comp=$6
snkRow=$7
srcRow=$8

# Summed ratio correlator time series
tmin=4
tstep=2
tmax=14

jkcfgs=$(( $cfgs - 1 ))
# The actual number of tseps per SR correlator
numT=$(( ($tmax-$tmin)/$tstep + 1 ))
# The full tseries of SR correlator
tseriesSR="${tmin}.${tstep}.${tmax}"
echo "SR data computed in range : $tseriesSR"

fitModel="summedRatio" # M + Cexp(-dE*T) + D*Texp(-dE*T)
tseriesSRFit="${tminfit}.${tstepfit}.${tmaxfit}" # fit summed ratio correlator according to: <min>.<step>.<max>
echo "Fitting SR data in range : $tseriesSRFit"

###########################
# Now do the jackknifing and xmbf parsing for the SR jackknife files
###########################
echo ""
echo "**********Working on SR correlator***********"
echo ""
# # Isolate base SR name by removing preceeding directory information
# cp $nptCorrsFile ${nptCorrsFile}.bak
# slashcount=`head -1 $nptCorrsFile | awk -F"/" '{printf NF"\n"}'`

echo "  Starting generation of jackknife samples/xmbf inputs for different SR correlators"
# pushd SR > /dev/null
~/scripts/jackknife/make_jackknife_samples $sRFile JACK $numT $cfgs
# popd > /dev/null


####################################################
# NOW MAKE XMBF INPUTS FOR EACH JACKKNIFE SAMPLE
####################################################
echo "`pwd`/$sRFile" >> `pwd`/thisCorr.list
# File containing full path/name to correlator
fileNpt=`pwd`/thisCorr.list

cd JACK/
for j in `seq 0 1 $jkcfgs`; do


    outPrefix=`echo $sRFile | sed -e 's@.dat@'_jack${j}'@'`
    echo $outPrefix

    echo "`pwd`/${outPrefix}.dat" > tmp.list

    # Transform this jackknife dat file into a form interpretable by xmbf
    ~/src/simulfts.py -c $jkcfgs -x $comp -T $tseriesSR -f tmp.list -n 2 \
    	-o $outPrefix -s 1 -p 2

    parsedNpt="${outPrefix}_simul2pt.${comp}_XMBF.dat"

    # Generate xml input for xmbf
    ~/src/xmbfini.py -d $parsedNpt -D $fitModel -b 1 -B 0 \
    	-t $tseriesSRFit >& SR.${comp}.jack${j}.ini.xml


    # Do the fit, catching output
    warn=""
    ${HOME}/scripts/xmbf_scripts/run_XMBF.csh . SR.${comp}.jack${j}.ini.xml res/ cov/ #>& jack${j}.log
    # tail -f jack${j}.log

    # warn1=`grep WARNING jack${j}.log`
    # warn2=`grep nan jack${j}.log`
    # warn=`cat $warn1 $warn2`


    # # Search log file for errors/warnings/nans
    # if [ $warn != "" ]; then
    # 	echo "Found something not finite or dangerous..."
    # 	echo "Noting jackknife sample and moving on..."
    # 	echo $j >> BAD.list
    # else # If no errors or warnings, dump the parsed xmbf input and input xml
    # 	rm *XMBF* *.xml
    # fi
    rm tmp.list

done

rm $fileNpt

exit 0
