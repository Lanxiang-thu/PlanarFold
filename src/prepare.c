#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dsfly.h"

//-----------------------------------------------------------------------------
// InitControl: set default value to control parameters
//-----------------------------------------------------------------------------
void InitAll(control *ctl, prmtop *top, energy *eng, state *sdt, ncfile *ncf){
    ctl->ntpr = 1000;       // Writing frequency of output diagnostics file (mdout)
    ctl->ntwx = 1000;       // Writing frequency of coordinate trajectory file (mdcrd)
    ctl->ntr = 0;        // Flag for restraining specified atoms in Cartesian
    // SMD related
    ctl->ismd = 0;       // option for Steered MD
    ctl->dumpfreq = 1000;
    ctl->stretch_force = 1.0;   // the stretching force for c.f.SMD, unit in pN
    ctl->stretch_spring = 0.01;  // the spring constant for c.v.SMD, unit in kcal/mol/A
    ctl->dist_init = 10.0;
    ctl->dist_final = 100.0;
    ctl->restraint_wt = 0.1;    // Force constant for ntr=1 use
    ctl->restraint_eq = 10.0;    // Equil value for ntr=1 use
    ctl->rcut = 11.5;
    //
    ctl->nscm = 1000;      // Frequency of removing translation and rotation of COM
    ctl->nstlim = 10000;     // Number of MD-steps to be performed
    ctl->dt = 0.004;      // Time steps (ps)
    // temperature control
    ctl->ntt = 1;        // Thermostat type, 0/1 for no/Berendsen Thermostat
    ctl->temp0 = NULL;   // Reference temperature for ntt>0
    ctl->tempi = NULL;   // Initial temperature
    ctl->duration = NULL;   // The duration of thermostat
    ctl->ntemp = 1;     // the section of different thermostating stage
    ctl->tautp = 0.005;   // Time constant for temperature fluctuations (in ps)
    //
    ctl->iCollide = 0;   // Time constant for temperature fluctuations (in ps)
    
    // energy
    eng->E_bond = 0.0;          // Energy due to bond-stretching terms
    eng->E_angl = 0.0;          // Energy due to angle-bending terms
    eng->E_dihe = 0.0;          // Energy due to dihedral terms
    eng->E_vdw = 0.0;           // Energy due to van der Waarl terms
    eng->E_stacking = 0.0;      // Energy due to base-stacking terms
    eng->E_basepair = 0.0;      // Energy due to base-pairing terms
    eng->E_restraint = 0.0;     // Energy due to restraint terms
    eng->Epot = 0.0;            // Total potential energy
    eng->Ekin = 0.0;            // Total kinetic energy
    eng->Etot = 0.0;            // Total energy of the system

    // STATES
    sdt->crd = NULL;        // The (current) coordinates of the atoms
    sdt->vel = NULL;        // The (current) velocities of the atoms
    sdt->acl = NULL;        // The (current) acceleration on atoms
    sdt->NotSearched = NULL;         // The Flag denoting whether this nucleotide is Base-Pair-saturate

    // TOPOLOGY
    top->pointers = NULL;      // POINTERS
    top->charges = NULL;     // CHARGE
    top->mass = NULL;       // MASS
    top->iac = NULL;           // ATOM_TYPE_INDEX
    top->ico = NULL;           // NONBONDED_PARM_INDEX
    // Force field parameters read from the input file
    top->stk_k = NULL;      // Stacking force constant (Turner's 36 parms)
    top->bpFLAG = NULL;     // The number of HB for each
    //
    top->anglist = NULL;       // ANGLES_WITHOUT_HYDROGEN
    top->inv_mass = NULL;   // The inverse atomic masses                                    // need to be supplemented before gradient()
    top->cn1 = NULL;        // LENNARD_JONES_ACOEF
    top->fcn1 = NULL;       // LENNARD_JONES_ACOEF used for force calculation(*12)
    top->cn2 = NULL;        // LENNARD_JONES_BCOEF
    top->fcn2 = NULL;       // LENNARD_JONES_BCOEF used for force calculation(*6)
    top->BPlist5 = NULL;
    top->BPlist3 = NULL;
    top->aBPlist5 = NULL;
    top->aBPlist3 = NULL;
    top->cmap_max_row = NULL;
    top->nBP = 0.0;
    top->nStem = 0.0;
    top->Stem = NULL;
    top->aStem = NULL;
    top->isTermi = NULL;
    top->BPlistTp = NULL;
    top->stklist = NULL;
    top->nSTk = 0.0;
    top->isStem = NULL;
    // restrainted MD
    top->ntrBplist = NULL;
    top->ntr_nBP = 0.0;
    top->dprBplist = NULL;
    top->dpr_nBP = 0.0;
    top->rcut = 0.0;
    // targeted MD
    top->tgtBplist = NULL;
    top->notTargeted = NULL;
    top->n_tgtBP = 0.0;
    top->orient = NULL;
    //SMD-related
    top->stretch_force = 0.0;
    top->dist_init = 0.0;
    top->dist_final = 0.0;

    // Restart file
    ncf->ncid_Rst = 0;
    ncf->varid_time_Rst = 0;     // The variable ID for 'time' of restart file
    ncf->varid_spatial_Rst = 0;  // The variable ID for 'spatial' of restart file
    ncf->varid_crd_Rst = 0;      // The variable ID for 'coordinates' of restart file
    ncf->varid_vel_Rst = 0;      // The variable ID for 'velocities' of restart file
    ncf->varid_cellSpatial_Rst = 0;  // The variable ID for 'cell spatial' of restart file
    ncf->varid_cellAngular_Rst = 0;  // The variable ID for 'cell_angular' of restart file
    ncf->varid_cellLengths_Rst = 0;  // The variable ID for 'cell_lengths' of restart file
    ncf->varid_cellAngles_Rst = 0;   // The variable ID for 'cell_angles' of restart file
}



//-----------------------------------------------------------------------------
// prepare: all the precomputation work during a single frame
//-----------------------------------------------------------------------------
void prepare(control *ctl, prmtop *top, energy *eng,
             state *sdt, timer *tim, FILE *fp, int nstep){
    int i, j, k, l, r;
    int n, natom, ntypes, nres;
    int idx, *ntyp, iparm, cont;
    int nex, *nexcld;
    int ii, jj, kk, n8;
    int nangl_H, nangl_X;
    int i3, ii2, i5, inew;
    int ia, ib, ic, id;
    int iaf, ibf, icf, idf;
    int c, d;
    int isfourier, bpdamp;
    int *dihe_fourier;
    int *buffer1, *dihedral;
    int *dihe_H, *dihe_X;
    int inbd;
    int ndihe, ndihe_H, ndihe_X;
    int *intra_nbd;
    int nparm;
    int ncell[3];
    int rboxes, rboxes2;
    int nboxes;
    int *boxvec;
    int min2;   //17
    int P;      //18
    int nwithq;     //5 charge round
    double tq;      //5 charge round
    double E_lvdw, deng_tail, invrc, invrc3, invrc9, ddeng;
    double rn[3];
    double maxbox;
    double *cn1, *cn2, *fcn1, *fcn2;
    double *sinabc, *cosabc;
    double *box;
    double *invm;
    double q, eele_corr;
    double PiVolume, PiBeta;   //18

    printf("data has entered prepare()\n");
    natom = top->pointers[0];
    sdt->natom = natom;
    ntypes = top->pointers[1];  // Number of distinct Lennard-Jones atom types
    nres = top->pointers[11];   // Number of residues
// use cubic box !                                               // delete this if using original box shape
    for (i = 3; i < 6; i++){
        sdt->box[i] = 90.0;
    }
    top->restraint_wt = ctl->restraint_wt;
    top->restraint_eq = ctl->restraint_eq;
//0 scaling all the potentials
    top->bond_k *= top->f_scale;
    top->epsilon *= top->f_scale;
    top->charge *= sqrt(top->f_scale);
    top->bphm_k *= top->f_scale;
    top->bphm_p *= top->f_scale;
    top->angle_k *= top->f_scale;
    top->angle_p *= top->f_scale;
    for (i = 0; i < 6*6; i++){
        top->stk_k[i] *= top->f_scale;
    }
    //top->stk_p *= top->f_scale;
//  verlet list
    sdt->verlet_list = (int **)malloc(natom*sizeof(int *));
    // the number of actual neighboring atoms
    sdt->verlet_natom =  (int *)malloc(natom * sizeof(int));
    for (i = 0; i < natom; i++){
        sdt->verlet_list[i] =  (int *)malloc(natom*sizeof(int));
    }
    

// Prepare the non-bonded list (exclude the residues without non-bonded interaction)
    int iFind;
    top->nbres = natom - top->nNBres;
    top->nblist =  (int *)malloc(top->nbres * sizeof(int));
    printf("top->inb:\n");
    for (i = 0; i < top->nNBres; i++){
        printf("%d\t", top->inb[i]);
    }
    n = 0;
    for (i = 0; i < natom; i++){
        iFind = 0;
        for (j = 0; j < top->nNBres; j++){
            if (i==top->inb[j]){
                iFind = 1;
                break;
            }
        }
        // not excluded
        if (iFind==0){
            top->nblist[n] = i;
            n += 1;
        }
    }
    printf("top->nblist:\n");
    for (i = 0; i < top->nbres; i++){
        printf("%d\t", top->nblist[i]);
    }

//1 perform dynamic allocation of some global arrays
    sdt->acl  = (double *)malloc(3*natom*sizeof(double));
    top->fcn1 = (double *)malloc((ntypes*(ntypes+1)/2)*sizeof(double));
    top->fcn2 = (double *)malloc((ntypes*(ntypes+1)/2)*sizeof(double));
    top->fcn2 = (double *)malloc((ntypes*(ntypes+1)/2)*sizeof(double));
    // 4 for atom[i,j,k] and angle type index
    // 4 means one bp-step has 4 angles
    top->inv_mass = (double *)malloc(natom*sizeof(double));
    //sdt->moved = (int *)malloc(natom*sizeof(int));
    //sdt->update = (int *)malloc(natom*sizeof(int));
    printf("\nprepare PART 1 is okay\n");
//3 transform the box angles form degree to radian
    for (i = 3; i < 6; i++){
        sdt->box[i] = sdt->box[i]* PI / 180.0;
    }
    printf("prepare PART 2 is okay\n");
//4 create inverse masses list
    double MassTot = 0.0;
    invm = top->inv_mass;
    for (i = 0; i < natom; i++){
        MassTot += top->mass[i];
        //printf("%lf\n", top->mass[i]);
        invm[i] = 1.0 / top->mass[i];
    }
    top->MassTot = MassTot;
    printf("prepare PART 3 is okay\n");

//14 zero out first derivatives (or 'forces')
    for (i = 0; i < 3*natom; i++){
        sdt->acl[i] = 0.0;
    }
    printf("prepare PART 4 is okay\n");

// Turner's parameter (GA-parameter scaling)
    for (i = 0; i < 6*6; i++){
        top->stk_k[i] *= top->stk_ks;
    }
    
// VDW
    top->ljA = 4 * top->epsilon * pow(top->sigma, 12);
    top->ljB = 4 * top->epsilon * pow(top->sigma, 6);
    top->fljA = 12.0 * top->ljA;
    top->fljB = 6.0 * top->ljB;
// BasePair
    /*
    for (i = 0; i < 10; i++){
        top->bpFLAG[i] = 1;
    }
    // index[AU]=2-1 and index[GC]=9-1
    top->bpFLAG[1] = 4;
    top->bpFLAG[4] = 4;
    top->bpFLAG[8] = 6;
    // G-A has been added!!! 2022-9-8
    top->bpFLAG[3] = 2;
    */
    // Dynamic programming use
    sdt->NotSearched = (int **)malloc(natom*sizeof(int *));
    sdt->cMap =   (double **)malloc(natom*sizeof(double *));
    top->cmap_max_row = (double *)malloc(natom*sizeof(double));
    sdt->triMap = (double **)malloc(natom*sizeof(double *));
    sdt->sumMap = (double **)malloc(natom*sizeof(double *));
    sdt->col_iblock = (int **)malloc(natom*sizeof(int *));
    sdt->row_iblock = (int **)malloc(natom*sizeof(int *));
    for (i = 0; i < natom; i++){
        sdt->NotSearched[i] =  (int *)malloc(natom*sizeof(int));
        sdt->cMap[i] =   (double *)malloc(natom*sizeof(double));
        sdt->triMap[i] = (double *)malloc(natom*sizeof(double));
        sdt->sumMap[i] = (double *)malloc(natom*sizeof(double));
        sdt->col_iblock[i] =  (int *)malloc(natom*sizeof(int));
        sdt->row_iblock[i] =  (int *)malloc(natom*sizeof(int));
    }
    // the block where resides the maximum cmap
    sdt->max_iblock = (int *)malloc(natom*sizeof(int));
    sdt->stlist = (int *)malloc(natom*sizeof(int));
    sdt->ndlist = (int *)malloc(natom*sizeof(int));
    // store the index+1 of the pairing residue, if 0 then unpairs
    // Dimension: 2 x natom
    top->BP_partners = (int **)malloc(2*sizeof(int *));
    for (i = 0; i < 2; i++){
        top->BP_partners[i] =  (int *)malloc(natom*sizeof(int));
    }
    top->nBP_partner = (int *)malloc(natom*sizeof(int));
    //
    top->aBPlist5 = (int *)malloc(natom*sizeof(int));
    top->aBPlist3 = (int *)malloc(natom*sizeof(int));
    top->BPlist5 = (int *)malloc(natom*sizeof(int));
    top->BPlist3 = (int *)malloc(natom*sizeof(int));    
    top->BPlistTp = (int *)malloc(natom*sizeof(int));
    // targeted MD    
    top->tgtBplist = (int *)malloc(2*natom*sizeof(int));
    top->notTargeted = (int *)malloc(natom*sizeof(int));
    top->n_tgtBP = 0;
    for (i = 0; i < natom; i++){
        top->notTargeted[i] = 1;
    }
    //
    top->nBP = 0;
    top->nStem = 0;
    top->aStem = (int *)malloc(natom*sizeof(int));
    top->Stem = (int *)malloc(natom*sizeof(int));
    top->isTermi = (int *)malloc(natom*sizeof(int));
    top->isStem = (int *)malloc(natom*sizeof(int));
    // For Pseudoknot-Exclusion use
    top->orient = (int *)malloc(natom*sizeof(int));

// Base Stacking
    top->anglist = (int *)malloc(2*natom*sizeof(int));
    top->stklist = (int *)malloc(4*4*natom*sizeof(int));
    top->nSTk = 0;

// random number (for water collison)
    sdt->rndcon = ctl->iCollide;
    // we always keep the random seed as 1.
    // in the future, we might adopt timer as random seed

//21 prepare for Thermostat
    // Determine the number of system degrees of freedom and particles (3N-3-Nc, I don't know why '-3' as well...)
    sdt->ndf = 2*natom;
    //sdt->Ttarget = ctl->tempi[0];
    printf("degree of freedom: ndf=%d\n", sdt->ndf);
    printf("prepare PART 5 is okay\n");


///===========================================================
    // c.f. SMD usage
    //sdt->spring_len = 

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
tim->t_Prepare += punchClock(tim->tt);

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//20 update the vel(-0.5) to vel(0) using acl(0) acquired by crd(0)
    const double dt = ctl->dt;  // sqrt(418.4) or just 418.4??              Only weird mdgx codes will add 418.4 here
    const double hdt = 0.5*dt;
    // convert velocity from half unit to complete unit (form kcal to Joule)
    i3 = 0;
    /*
    for (i = 0; i < natom; i++){
        //printf("%8.4lf  %8.4lf \n", sdt->vel[i3], sdt->vel[i3+1]);
        sdt->vel[i3] *= 20.4548282808730;   // 20.4548282808730 = srqt(/418.4)
        sdt->vel[i3+1] *= 20.4548282808730;
        i3 += 3;
    }
    */
tim->t_Verlet += punchClock(tim->tt);
    printf("prepare PART 6 is okay\n");

// "mdgx regime"
// Generate the verlet list (every nVerlet steps)
    
    int ia3, ib3;
    double dx, dy, r2;
    double vlcut2 = ctl->verlet_cut * ctl->verlet_cut;
    int **verlet_list = sdt->verlet_list;
    int *verlet_natom = sdt->verlet_natom;
    // enlarge or shrink the verlet list according to the instantaneous Temperature
    // (the default temperature is 300K: if Temp(t)=900K, vlcut = 3*vlcut)
    if (sdt->T > 300.0){ 
        vlcut2 *= sqrt(sdt->T / 300.0);
    }
    //
    for (i = 0; i < natom-2; i++){
        n = 0;
        ia3 = i*3;
        verlet_natom[i] = 0;
        // j=i+2 as we exclude two adjacent nucleotides
        // G--A--U--C--'G'--A--C--G--U
        // o--o--o--x--'*'--x--o--o--o
        for (j = i+2; j < natom; j++){
            ib3 = j*3;
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > vlcut2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            else{
                // store the neighboring residues into verlet buffer "verlet_list"
                n = verlet_natom[i];
                verlet_list[i][n] = j;    // store the index of coordinates (x3)
                verlet_natom[i] += 1;
            }
        }
    }
    
    //  get the potential energy and update atomic forces at full step, acl(0)
    // zero out first derivatives (or 'forces')
    for (i = 0; i < 3*natom; i++){
        sdt->acl[i] = 0.0;
    }
    gradient(eng, top, sdt, tim, ctl, fp, nstep);
    printf("prepare PART 6-2 is okay\n");
    // use Newton second law to get the next accelerations;
    // find the full-step velocities using the Verlet recursion
    //printf("dt=%lf\n", dt);
    i3 = 0;
    //printf("Full-step velocities\n");
    double invmass;
    for (i = 0; i < natom; i++){
        invmass = top->inv_mass[i] * -418.4;           // 418.4 will convert kcal to J
        for (j = 0; j < 2; j++){
            // now we get the real 'acceleration' from former cumputed
            // first derivatives in gradient()
            sdt->acl[i3+j] *= invmass;
        }
        //printf("%d, %lf, %lf\t", i, sdt->vel[i3], sdt->vel[i3+1]);
        i3 += 3;
    }
    //------test ------------------------------------------------------------------------
    /*
    i3 = 0;
    for (i = 0; i < natom; i++){
        printf("%d, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf\n", i, sdt->crd[i3], sdt->crd[i3+1], sdt->vel[i3], sdt->vel[i3+1], sdt->acl[i3], sdt->acl[i3+1]);
        i3 += 3;
    }
    */
    //------test ------------------------------------------------------------------------
    /*
    i3 = 0;
    double ivm;
    for (i = 0; i < natom; i++){
        //printf("crd:%8.4lf  %8.4lf \t", sdt->crd[i3], sdt->crd[i3+1]);
        //printf("vel:%8.4lf  %8.4lf \n", sdt->vel[i3], sdt->vel[i3+1]);
        ivm = invm[i] * -418.4;           // 418.4 will convert kcal to J
        // now we get the real 'acceleration' from former cumputed
        // first derivatives in gradient()
        sdt->acl[i3] *= ivm;
        sdt->acl[i3+1] *= ivm;
        //sdt->acl[i3+2] *= ivm;
        // get the full-step velocities
        sdt->vel[i3] += sdt->acl[i3]*hdt;
        sdt->vel[i3+1] += sdt->acl[i3+1]*hdt;
        //sdt->vel[i3+2] += sdt->acl[i3+2]*hdt;
        i3 += 3;
    }
    */
    
    
tim->t_Verlet += punchClock(tim->tt);
    printf("prepare PART 7 is okay\n");
    


//22 Thermostat
    // Compute kinetic energy and instantaneous temperature after the force computation  
    /*
    double velSum;
    double tscale;
    double lambda;
    i3 = 0;
    // kinetic energy: Ek
    eng->Ekin = 0.0;    // zero out kinetic energy
    for (i = 0; i < natom; i++){
        velSum = sdt->vel[i3] * sdt->vel[i3];
        velSum += sdt->vel[i3+1] * sdt->vel[i3+1];
        //velSum += sdt->vel[i3+2] * sdt->vel[i3+2];
        eng->Ekin +=  top->mass[i] * velSum;
        i3 += 3;
    }
    eng->Ekin /= 836.8;     // eng->Ekin = 0.5*eng->Ekin /418.4
    eng->Etot = eng->Ekin + eng->Epot;
    // instantaneous temperature: T
    sdt->T = 2.0 * eng->Ekin / (sdt->ndf * GASCNST);
    
    // Berendsen thermostat
    if (ctl->ntt == 1){  
        // Temperature was computed after the force computation;
        // velocities have been updated since the temperature   
        // computation but that is part of the algorithm.
        if (sdt->T < 0.5*sdt->Ttarget) {
            tscale = 1.0;
        }
        // Too hot!  Let's call it twice the target temperature
        else if (sdt->T > 2.0*sdt->Ttarget) {
            tscale = -0.5;
        }
        // Ahh, just right...
        else {
            tscale = (sdt->Ttarget / sdt->T) - 1.0;
        }
        // get the scaling factor 'lambda'
        lambda = sqrt(1.0 + tscale * dt / ctl->tautp);
        // velocity scaling (no virtual sites)
        for (i = 0; i < 3*natom; i++){
            sdt->vel[i] *= lambda;
        }
    }
    */
    
    
tim->t_Thermostat += punchClock(tim->tt);
    printf("//-----------------------------------------------------------------------------\n");
    printf("Frame:   %12d   Temperature:    %12.4lf   Restraint:  %12.4lf\n", 0, sdt->T, eng->E_restraint);
    printf("Bond:    %12.4lf   Angle:          %12.4lf   Stacking:   %12.4lf\n", eng->E_bond, eng->E_angl, eng->E_stacking);
    printf("Vdw:     %12.4lf   Elec:           %12.4lf   BasePair:   %12.4lf\n", eng->E_vdw, eng->E_elec, eng->E_basepair);
    printf("Ektot:   %12.4lf   Eptot:          %12.4lf   Etot:       %12.4lf\n", eng->Ekin, eng->Epot, eng->Etot);
    //printf("!!!Ekin=%lf, Etot=%lf\n", eng->Ekin, eng->Etot);
    i3 = 0;
    /*
    for (i = 0; i < natom; i++){
        //printf("crd:%8.4lf  %8.4lf \t", sdt->crd[i3], sdt->crd[i3+1]);
        //printf("vel:%8.4lf  %8.4lf \n", sdt->vel[i3], sdt->vel[i3+1]);
        for (j = 0; j < 2; j++){
            if (isnan(sdt->vel[i3+j])){
                printf("atom-%d has velociy(%d): %lf\n", i, j, sdt->vel[i3+j]);
                exit(0);
            }
        }
        i3 += 3;
    }
    */
    printf("prepare PART 22 is okay\n");
}


