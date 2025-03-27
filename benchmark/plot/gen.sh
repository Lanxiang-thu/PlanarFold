#!/bin/bash

# generate topology and coordinate files
generate_PlanarFold_Inputs -s ssRNA.seq -c ssRNA.inpcrd -p ssRNA.prmtop\
			   -ntr restraint_ntr.in -dpr restraint_dpr.in \
			   --iCirc 1 --InitV 0 --vscale 0.0



