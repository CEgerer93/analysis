#!/usr/bin/bash

## Use this fella to form summed ratios for various mod(pz) & mod(zsep)

if [ $# -ne 5 ]; then
    echo "Usage: $0 <pz> <zsep> <gamma> <SR fit tmin> <SR fit tmax>"
    exit 1
fi

phaseStub=unphased
operator=''

op=''
oprow=0

tmin3pt=$4
tstep=2
tmax3pt=$5

if [ $3 -eq 8 ]; then
    op=b_b0xDA__J0_A1pP
    oprow=1
elif [ $3 -eq 11 ]; then
    op=a_a1xDA__J1_T1pM
    oprow=2
fi

# Get correct 2pt src tslice
tslice2pt=0
if [ $1 -eq 0 ] || [ $1 -eq 1 ]; then
    tslice2pt=-2
fi

# Get correct num tslices of 2pt
tmax2pt=15
if [ $1 -eq 0 ] || [ $1 -eq 1 ] || [ $1 -eq 4 ]; then
    tmax2pt=20
fi

# Modulus of momentum
modmom=`echo "sqrt($1*$1)" | bc`

for p in $1; do

    # if [ $p -gt 3 ]; then
    # 	phaseStub=phased
    # fi    

    for z in $2; do

	# Remove old sum.pz.xml file
	rm sum.${p}${z}.xml
	
	disp_list=""
	for n in `seq 1 1 $z`; do
	    disp_list="${disp_list}3 "
	done

	
	# For standard 3pt/2pt ratios w/ pz != 0
	sed -e 's@PZ@'"${p}"'@g' -e 's@MODMOM@'"${modmom}"'@g' -e 's@DZ@'"${disp_list}"'@g' \
	    -e 's@TMAX@'"${tmax2pt}"'@g' -e 's@TSLICE2PT@'"${tslice2pt}"'@g' \
	    -e 's@ROWINS@'"${oprow}"'@g' -e 's@INS@'"${op}"'@g' \
	    -e 's@3PTTmin@'"${tmin3pt}"'@g' -e 's@3PTTmax@'"${tmax3pt}"'@g' \
	    sum-matelem.ini.xml > sum.${p}${z}.xml

	    # -e 's@XPHASEX@'"${phaseStub}"'@g' -e 's@XOPX@'"${operator}"'@g' \
	    # sum-matelem.ini.xml > sum.${p}${z}.xml

	./pITD-matelem sum.${p}${z}.xml
	
	
    done
done
