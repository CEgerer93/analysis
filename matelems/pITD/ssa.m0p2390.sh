#!/usr/bin/bash

## Use this fella to form summed ratios for various mod(pz) & mod(zsep)

if [ $# -ne 5 ]; then
    echo "Usage: $0 <pz> <zsep> <gamma> <SR fit tmin> <SR fit tmax>"
    exit 1
fi

EXE=/w/scifs17/JLabLQCD/cegerer/pPDF/pITD-matelem
XML=/work/JLabLQCD/cegerer/pPDF/sum-matelem.m0p2390.ini.xml

PHASEDIR=unphased
PHASESTUB=$PHASEDIR
if [ $1 -gt 3 ] || [ $1 -lt -3 ]; then
    PHASEDIR=phased
fi
if [ $1 -gt 3 ]; then
    PHASESTUB=${PHASEDIR}/d001_2.00
elif [ $1 -lt -3 ]; then
    PHASESTUB=${PHASEDIR}/d001_-2.00
fi

operator=''

op=''
oprow=0
gamma=0
PROJ=0

##############################################################################
# Set the 2pt fitting window - sloppy for now!
tmin2ptFit=2; tmax2ptFit=0
if [ $1 -eq 0 ]; then tmax2ptFit=18; fi
if [ $1 -eq 1 ] || [ $1 -eq -1 ]; then tmax2ptFit=15; fi
if [ $1 -eq 2 ]; then tmax2ptFit=12; fi
if [ $1 -eq -2 ]; then tmax2ptFit=10; fi
if [ $1 -eq 3 ] || [ $1 -eq -3 ]; then tmax2ptFit=11; fi
if [ $1 -eq 4 ] || [ $1 -eq -4 ]; then tmax2ptFit=11; fi
if [ $1 -eq 5 ]; then tmax2ptFit=11; fi
if [ $1 -eq -5 ]; then tmax2ptFit=10; fi
if [ $1 -eq 6 ]; then tmax2ptFit=8; fi
if [ $1 -eq -6 ]; then tmax2ptFit=8; fi

# Set the 3pt fitting window
tmin3ptFit=$4
tstep=2
tmax3ptFit=$5
#------------------------------------------------------------------------------

if [ $3 -eq 8 ]; then
    op=b_b0xDA__J0_A1pP
    oprow=1
    gamma=8
    PROJ=1
elif [ $3 -eq 11 ]; then
    op=a_a1xDA__J1_T1pM
    oprow=2
    gamma=11
    PROJ=2
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
	    -e 's@3PTTmin@'"${tmin3ptFit}"'@g' -e 's@3PTTmax@'"${tmax3ptFit}"'@g' \
	    -e 's@2PTTmin@'"${tmin2ptFit}"'@g' -e 's@2PTTmax@'"${tmax2ptFit}"'@g' \
	    -e 's@PHASEDIR@'"${PHASEDIR}"'@g' -e 's@PHASESTUB@'"${PHASESTUB}"'@g' \
	    -e 's@GAMMA@'"${gamma}"'@' -e 's@PROJ@'"${PROJ}"'@' \
	    $XML > sum.${p}${z}.xml

	    # -e 's@XPHASEX@'"${phaseStub}"'@g' -e 's@XOPX@'"${operator}"'@g' \
	    # sum-matelem.ini.xml > sum.${p}${z}.xml

	$EXE sum.${p}${z}.xml
	
	
    done
done
