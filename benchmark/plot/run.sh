#!/bin/bash

# Step-1: Filter the secondary structure to keep only apical stems
extract_Apical_loops --seq ssRNA.seq  --seq_apical ssRNA_apical.seq

# Step-2: Generate restraint and confinement files for apical stems
generate_PlanarFold_Inputs -s ssRNA_apical.seq  -ntr restraint_ntr_apical.in -dpr restraint_dpr_apical.in

# Step-3: Pre-equilibrium to restrain the apical stems to avoid trans-isomer
PlanarFold  -i controls_prep.in  -f fparm.in -p ssRNA.prmtop  -c ssRNA.inpcrd \
  	    -o run_prep.out  -x run_prep.nc  -r run_prep.rst \
	    -ntr restraint_ntr_apical.in -dpr restraint_dpr_apical.in


# Step-4: Production run with all restraints applied
PlanarFold  -i controls.in  -f fparm.in -p ssRNA.prmtop  -c run_prep.rst \
  	    -o run.out  -x run.nc  -r run.rst \
	    -ntr restraint_ntr.in -dpr restraint_dpr.in





