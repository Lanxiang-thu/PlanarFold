#!/bin/bash

# plot
plot_Frames -o run.out -p ssRNA.prmtop -x run.nc \
	    --frame0 25 --frame1 26 --step 25\
	    --needRotate 1 -s 5 --needEnergy 1



