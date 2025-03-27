#!/bin/bash

# start simulation using ssRNA.inpcrd
#PlanarFold  -i controls.in  -f fparm.in -p ssRNA.prmtop  -c ssRNA.inpcrd \
#  	      -o run001.out  -x run001.nc  -r run001.rst \
#	      -ntr restraint_ntr.in -dpr restraint_dpr.in
#cat dist.RST > dist001.RST

# extend simulation using restart file
istart=2
iend=10
I=$istart
while [ $I -le $iend ]
do
  suffix=`printf %03d $I`
  suffix2=`printf %03d $((I-1))`
  # production run
  PlanarFold  -i controls.in  -f fparm.in -p ssRNA.prmtop  -c run$suffix2.rst \
  	      -o run$suffix.out  -x run$suffix.nc  -r run$suffix.rst \
	      -ntr restraint_ntr.in -dpr restraint_dpr.in
  cat dist.RST > dist$suffix.RST  
  I=$((I+1))
done


