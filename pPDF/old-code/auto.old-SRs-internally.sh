#!/usr/bin/bash


dataTMin=4
dataTMax=14

# Actual fit tseries
tminfit=4
tstepfit=2
tmaxfit=14

# 0-based!
excludedTsFromFit=-1
# excludedTsFromFit=0
# excludedTsFromFit=4,5
# excludedTsFromFit=5


### File prefix
# pref="tS"
pref="tm60"

### Norm type
normType="rest-norm"
# normType="momx-norm"
# normType="p001-norm"


# suffix=summedRatio_p300
# suffix=summedRatio_p001


cfgs=349
jkcfgs=$(( $cfgs - 1 ))

srcrow=1
snkrow=1
# insertion=pion_pion_2xDA__J0_A1mM
insertion=b_b0xDA__J0_A1pP
# insertion=a_a1xDA__J1_T1pM
insrow=1


# Here for doing the fits
for comp in 1 2; do
   for f in ${pref}*.dat; do
    
	# ./ftsSR.sh $f <min tsep> <tsep step> <max tsep> <cfgs> <complexity (Re-1) or (Im-2)> <snk row> <src row>
	./ftsSR.sh $f $tminfit $tstepfit $tmaxfit $cfgs $comp 1 1

	resfile=`echo $f | sed -e 's@.dat@'".comp_${comp}.RES.dat"'@'`
	covfile=`echo $f | sed -e 's@.dat@'".comp_${comp}.COV.dat"'@'`

	# Now concatenate the results
	rm RES/$resfile COV/$covfile
	~/structure/rescat.csh $jkcfgs RES/$resfile COV/$covfile

	
    done
done


# # Here for gathering the matelems
# for p in 1 2 3; do
#     suffix=summedRatio_p${p}00
# # for p in `seq 1 1 6`; do
# # for p in `seq 3 1 3`; do
#     for comp in 1 2; do
# 	# for z in 8; do
# 	for z in `seq 0 1 8`; do
# 	# for z in `seq 9 1 16`; do
	    
# 	    disp=""
# 	    if [ $z -eq 0 ]; then
# 		disp="."
# 	    elif [ $z -lt 0 ]; then
# 		disp=","
# 		tmpZ=$(( $z * -1 ))
# 		for n in `seq 1 1 $tmpZ`; do
# 		    disp="${disp}-3"
# 		done
# 		disp="${disp}."
# 	    elif [ $z -gt 0 ]; then
# 		disp=","
# 		for n in `seq 1 1 $z`; do
# 		    disp="${disp}3"
# 		done
# 		disp="${disp}."
# 	    fi
	    
	    
	    
# 	    ./SRfits.py -c $cfgs -x $comp -t 6 -z $z -p $p \
# 	    	-n SRs/${normType}/${insertion}/${pref},fI1Y3i1,r${snkrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.tm3,fI2Y0i0,r${insrow},000,${insertion}__000${disp}t0,fI1Y3i1,r${srcrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.${suffix}.dat \
# 	    	-r RES/${pref},fI1Y3i1,r${snkrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.tm3,fI2Y0i0,r${insrow},000,${insertion}__000${disp}t0,fI1Y3i1,r${srcrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.${suffix}.comp_${comp}.RES.dat \
# 	    	-v COV/${pref},fI1Y3i1,r${snkrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.tm3,fI2Y0i0,r${insrow},000,${insertion}__000${disp}t0,fI1Y3i1,r${srcrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.${suffix}.comp_${comp}.COV.dat \
# 	    	-I $insertion -w $dataTMin -y $dataTMax -E $excludedTsFromFit

# 	    ##########
# 	    # For without summedRatio_pXYZ suffix
# 	    ##########
# 	    # ./SRfits.py -c $cfgs -x $comp -t 6 -z $z -p $p \
#             #     -n SRs/${normType}/${insertion}/${pref},fI1Y3i1,r${snkrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.tm3,fI2Y0i0,r${insrow},000,${insertion}__000${disp}t0,fI1Y3i1,r${srcrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.summedRatio.dat \
#             #     -r RES/${pref},fI1Y3i1,r${snkrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.tm3,fI2Y0i0,r${insrow},000,${insertion}__000${disp}t0,fI1Y3i1,r${srcrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.summedRatio.comp_${comp}.RES.dat \
#             #     -v COV/${pref},fI1Y3i1,r${snkrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.tm3,fI2Y0i0,r${insrow},000,${insertion}__000${disp}t0,fI1Y3i1,r${srcrow},00${p},NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1__${p}00.summedRatio.comp_${comp}.COV.dat \
#             #     -I $insertion -w $dataTMin -y $dataTMax -E $excludedTsFromFit
	    
# 	done
#     done
# done
