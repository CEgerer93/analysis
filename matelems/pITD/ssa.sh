#!/usr/bin/bash

## Use this fella to form summed ratios for various mod(pz) & mod(zsep)

phaseStub=unphased
operator=''

# for p in `seq 4 1 6`; do
# for p in 0 1 2 3; do
for p in 0; do

    if [ $p -gt 3 ]; then
	phaseStub=phased
    fi

    if [ $p -eq 0 ]; then
	operator=NucleonMG1g1MxD0J0S_J1o2_G1g1
    elif [ $p -ne 0 ]; then
	operator=NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1
    fi
    

    for z in 1; do
    # for z in `seq 0 1 8`; do
    # for z in `seq 0 1 16`; do
    # for z in `seq 9 1 16`; do

	# Remove old sum.pz.xml file
	rm sum.${p}${z}.xml
	
	disp_list=""
	for n in `seq 1 1 $z`; do
	    disp_list="${disp_list}3 "
	done
	
	# sed -e 's@PZ@'"${p}"'@g' -e 's@DZ@'"${disp_list}"'@g' sum.arb-norm.tmp.ini.xml > sum.${p}${z}.xml
	# sed -e 's@PZ@'"${p}"'@g' -e 's@DZ@'"${disp_list}"'@g' sum.rest-norm.tmp.ini.xml > sum.${p}${z}.xml
	# sed -e 's@PZ@'"${p}"'@g' -e 's@DZ@'"${disp_list}"'@g' sum.pzneq0-norm.tmp.ini.xml > sum.${p}${z}.xml

	
	# For standard 3pt/2pt ratios w/ pz != 0
	sed -e 's@PZ@'"${p}"'@g' -e 's@DZ@'"${disp_list}"'@g' \
	    -e 's@XPHASEX@'"${phaseStub}"'@g' -e 's@XOPX@'"${operator}"'@g' \
	    sum-matelem.ini.xml > sum.${p}${z}.xml

	# # For standard 3pt/2pt ratios w/ px != 0
	# sed -e 's@PX@'"${p}"'@g' -e 's@DZ@'"${disp_list}"'@g' \
	#     -e 's@XPHASEX@'"${phaseStub}"'@g' sum-matelem_xmom.ini.xml > sum.${p}${z}.xml
	
	
	/home/cegerer/analysis/pPDF/summation-matelem sum.${p}${z}.xml
	# ./summation-matelem sum.${p}${z}.xml

	# mv *RE_fit* /work/JLabLQCD/cegerer/pPDF/Summed-Matelems/b_b0xDA__J0_A1pP/RES/
	# mv *IM_fit* /work/JLabLQCD/cegerer/pPDF/Summed-Matelems/b_b0xDA__J0_A1pP/RES/

	# mv tm60,fI1Y3i1,r1,00* Summed-Matelems/b_b0xDA__J0_A1pP/

	# rm *RE_fit*
	# rm *IM_fit*
	#rm tm60,fI1Y3i1,r1,00*
	
	
    done
done
