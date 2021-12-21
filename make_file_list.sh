#!/bin/bash 

if [ $# -ne 4 ]; then
    echo "Usage: $0 <pz> <chroma gamma> <phase (0 -or- 1)> <ensem>"
    exit 1
fi

pz=$1
gamma=$2
phase=$3
ensem=$4

PREFIX=''
if [ $phase -eq 0 ]
then
    PREFIX=/cache/lqcdpdf/isoClover/dist_ppdf/${ensem}/unphased/t0_avg
elif [ $phase -eq 1 ]
then
    PREFIX=/cache/lqcdpdf/isoClover/dist_ppdf/${ensem}/phased/t0_avg
fi


# Acceptable operator names
rest=NucleonMG1g1MxD0J0S_J1o2_G1g1
pboost=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pz1=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pzm1=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pz2=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pzm2=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pz3=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pzm3=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pz4=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pzm4=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pz5=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pzm5=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pz6=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# pzm6=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
# Set the operators
Op=''


if [ $pz -eq 0 ]; then Op=$rest; fi
if [ $pz -ne 0 ]; then Op=$pboost; fi
# if [ $pz == -1 ]; then Op=$pzm1; fi
# if [ $pz == 2 ]; then Op=$pz2; fi
# if [ $pz == -2 ]; then Op=$pzm2; fi
# if [ $pz == 3 ]; then Op=$pz3; fi
# if [ $pz == -3 ]; then Op=$pzm3; fi
# if [ $pz == 4 ]; then Op=$pz4; fi
# if [ $pz == -4 ]; then Op=$pzm4; fi
# if [ $pz == 5 ]; then Op=$pz5; fi
# if [ $pz == -5 ]; then Op=$pzm5; fi
# if [ $pz == 6 ]; then Op=$pz6; fi
# if [ $pz == -6 ]; then Op=$pzm6; fi


# Modulus of in/out-going momenta
modPx=0
modPy=0
modPz=`echo "sqrt($pz*$pz)" | bc `

# Determine the star(p)
for pi in $modPx $modPy $modPz; do
    echo $pi >> poo.list
done
pstar=`sort -rn poo.list | xargs | awk -F' ' '{printf $1$2$3"\n"}'`
rm poo.list



# Determine modulus of momentum transfer
qx=0
qy=0
qz=0
modQX=0
modQY=0
modQZ=0
qstar='000'




for z in `seq -8 1 8`; do

    ztag=''

    if [ $z -lt 0 ]; then
	for zi in `seq -1 -1 $z`; do
	    ztag="${ztag}-3"
	done
	ztag=",${ztag}"
    elif [ $z -gt 0 ]; then
	for zi in `seq 1 1 $z`; do
	    ztag="${ztag}3"
	done
	ztag=",${ztag}"
    else
	ztag=''
    fi

    # echo $ztag

    for r1 in 1 2; do
    	for r2 in 1 2; do
    	    for t in `seq 4 2 14`; do

		if [ $gamma -eq 8 ]; then
		    echo "${PREFIX}/tsnk_${t}/momXYZ.0.0.${pz}/ENS/t${t},fI1Y3i1,r${r2},00${pz},${Op}__${pstar}.tm3,fI2Y0i0,r1,000,b_b0xDA__J0_A1pP__${qstar}${ztag}.t0,fI1Y3i1,r${r1},00${pz},${Op}__${pstar}.dat" >> lists/b_b0xDA__J0_A1pP/r1/nuc3pt_p00${pz}_r${r2}${r1}_g8_z${z}.list
		fi

		if [ $gamma -eq 11 ]; then
    		    echo "${PREFIX}/tsnk_${t}/momXYZ.0.0.${pz}/ENS/t${t},fI1Y3i1,r${r2},00${pz},${Op}__${pstar}.tm3,fI2Y0i0,r2,000,a_a1xDA__J1_T1pM__${qstar}${ztag}.t0,fI1Y3i1,r${r1},00${pz},${Op}__${pstar}.dat" >> lists/a_a1xDA__J1_T1pM/r2/nuc3pt_p00${pz}_r${r2}${r1}_g11_z${z}.list
		fi

		if [ $gamma -eq 7 ]; then
		    echo "${PREFIX}/tsnk_${t}/momXYZ.0.0.${pz}/ENS/t${t},fI1Y3i1,r${r2},00${pz},${Op}__${pstar}.tm3,fI2Y0i0,r1,000,pion_pion_2xDA__J0_A1mM__${qstar}${ztag}.t0,fI1Y3i1,r${r1},00${pz},${Op}__${pstar}.dat" >> lists/pion_pion_2xDA__J0_A1mM/nuc3pt_p00${pz}_r${r2}${r1}_g7_z${z}.list
		fi

		if [ $gamma -eq 13 ] || [ $gamma -eq 14 ]; then
		    echo "${PREFIX}/tsnk_${t}/momXYZ.0.0.${pz}/ENS/t${t},fI1Y3i1,r${r2},00${pz},${Op}__${pstar}.tm3,fI2Y0i0,r1,000,a_a1xDA__J1_T1pM__${qstar}${ztag}.t0,fI1Y3i1,r${r1},00${pz},${Op}__${pstar}.dat" "${PREFIX}/tsnk${t}/momXYZ.0.0.${pz}/ENS/t${t},fI1Y3i1,r${r2},00${pz},${Op}__${pstar}.tm3,fI2Y0i0,r3,000,a_a1xDA__J1_T1pM__${qstar}${ztag}.t0,fI1Y3i1,r${r1},00${pz},${Op}__${pstar}.dat" >> lists/a_a1xDA__J1_T1pM/r1r3/nuc3pt_p00${pz}_r${r2}${r1}_g13andg14_z${z}.list
		fi

    	    done # t
    	done # r2
    done # r1

done # z

