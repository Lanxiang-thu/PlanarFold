#!/bin/bash

istart=10
iend=73

# Step1. prepare the initial unfolded conformation (skip this step if conformation is already provided)
cut -c1-$istart reverse_switch.seq > transcript.seq
generate_PlanarFold_Inputs -s transcript.seq -c ssRNA`printf %03d $istart`.inpcrd

# Step2. start the simulation
I=$istart
while [ $I -le $iend ]
do
  suffix=`printf %03d $I`
  suffix2=`printf %03d $((I+1))`

  # Generate the sequence file for this transcript
  cut -c1-$I reverse_switch.seq > transcript.seq

  # Generate topology file of this transcript using sequence file, but do not generate coordinate file
  generate_PlanarFold_Inputs -s transcript.seq -p ssRNA$suffix.prmtop -ntr restraint_ntr.in -dpr restraint_dpr.in

  # Production run
  PlanarFold  -i controls.in  -f fparm.in -p ssRNA$suffix.prmtop  -c ssRNA$suffix.inpcrd \
  	      -o run$suffix.out  -x run$suffix.nc  -r run$suffix.rst \
	      -ntr restraint_ntr.in -dpr restraint_dpr.in

  # Elongate the transcript (only adding 1-nt is allowed) using restart file, generate extended coordinate file
  elongate_Transcript -r run$suffix.rst -c ssRNA$suffix2.inpcrd

  I=$((I+1))
done

