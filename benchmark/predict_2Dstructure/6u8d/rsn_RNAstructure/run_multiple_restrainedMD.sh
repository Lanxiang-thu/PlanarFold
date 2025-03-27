#!/bin/bash

istart=1
iend=207

I=$istart
while [ $I -le $iend ]
do
  suffix=`printf %03d $I`
  # pre-equilibrium to release long-range restraints tension
  PlanarFold  -i controls_prep.in  -f fparm.in -p ssRNA.prmtop  -c ssRNA.inpcrd \
  	      -o run.out  -x run.nc  -r run_prep_$suffix.rst \
	      -ntr ./restraints/restraint_ntr_$suffix.in \
	      -dpr ./restraints/restraint_dpr_$suffix.in
  # production run
  PlanarFold  -i controls.in  -f fparm.in -p ssRNA.prmtop  -c run_prep_$suffix.rst \
  	      -o run$suffix.out  -x run$suffix.nc  -r run$suffix.rst \
	      -ntr ./restraints/restraint_ntr_$suffix.in \
	      -dpr ./restraints/restraint_dpr_$suffix.in
  I=$((I+1))
done




