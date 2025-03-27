#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <boost/python.hpp>
#include <numpy/ndarrayobject.h>
#include <time.h>
#include <iostream>
#include <typeinfo>     // test for variable type
#include "dsfly.h"

using namespace boost::python;
using namespace std;

void run_md(list toplist, list rstlist, list ctrlist, list prmlist, list ntrlist, list dprlist,
		    const std::string& fnOut, const std::string& fnTrj, const std::string& fnRst){
    int i, j, k, n, nstep, nframe;
    control ctl;
    prmtop top;
    energy eng;
    state sdt;
    ncfile ncf;
    timer tim;


// Initialize control parameters
    InitAll(&ctl, &top, &eng, &sdt, &ncf);
// Initialize timers
    InitTimer(&tim);

//-------------------------------------------------------------------------------------------------------------
// initialize force field parameters
    printf("Start: extract prmtop\n");
    // pointer
    object obj_pointers =    extract<object>(toplist[0]);
    top.pointers     =    (int *)(PyArray_DATA(obj_pointers.ptr()));
    // charge (in case of different charged residues)
    object obj_charges =    extract<object>(toplist[1]);
    top.charges     =    (double *)(PyArray_DATA(obj_charges.ptr()));
    // mass
    object obj_mass   =    extract<object>(toplist[2]);
    top.mass       = (double *)(PyArray_DATA(obj_mass.ptr()));
    // iac (ATOM_TYPE_INDEX)
    object obj_iac     =    extract<object>(toplist[3]);
    top.iac         =    (int *)(PyArray_DATA(obj_iac.ptr()));
    // ico (NONBONDED_PARM_INDEX)
    object obj_ico     =    extract<object>(toplist[4]);
    top.ico         =    (int *)(PyArray_DATA(obj_ico.ptr()));
    // inb (EXCLUDED_ATOMS_LIST)
    object obj_inb     =    extract<object>(toplist[5]);
    top.inb         =    (int *)(PyArray_DATA(obj_inb.ptr()));
    // nNBres
    top.nNBres         =    extract<int>(toplist[6]);

	printf("Done: extract prmtop\n");

//-----------------------------------------------------------------------------------------------------------
// initialize coordinates, velocities, and cell dimensions
    printf("Start: extract rst\n");
    // cell
    for (i = 0; i < 6; i++){
        sdt.box[i]   = extract<double>(rstlist[0][i]);
    }
    // crd
    object obj_crd     =    extract<object>(rstlist[1]);
    sdt.crd        = (double *)(PyArray_DATA(obj_crd.ptr()));
    // vel
    object obj_vel     =    extract<object>(rstlist[2]);
    sdt.vel        = (double *)(PyArray_DATA(obj_vel.ptr()));
	printf("Done: extract rst\n");

    // test for variable type length 
    cout<<"sizeof(sdt.crd[0]):"<<sizeof(sdt.crd[0])<<endl;
    cout<<"sizeof(sdt.vel[0]):"<<sizeof(sdt.vel[0])<<endl;
    //

//-----------------------------------------------------------------------------------------------------------
// initialize control variables
    printf("Start: extract control\n");
    ctl.ntpr            = extract<int>(ctrlist[0]);
    ctl.ntwx            = extract<int>(ctrlist[1]);
    ctl.ntr             = extract<int>(ctrlist[2]);
    //
    ctl.ismd            = extract<int>(ctrlist[3]);
    ctl.dumpfreq        = extract<int>(ctrlist[4]);
    ctl.stretch_force   = extract<double>(ctrlist[5]);
    ctl.stretch_spring  = extract<double>(ctrlist[6]);
    ctl.dist_init       = extract<double>(ctrlist[7]);
    ctl.dist_final      = extract<double>(ctrlist[8]);
    ctl.restraint_wt    = extract<double>(ctrlist[9]);
    ctl.restraint_eq    = extract<double>(ctrlist[10]);
    ctl.rcut            = extract<double>(ctrlist[11]);
    //
    ctl.nscm            = extract<int>(ctrlist[12]);
    ctl.nstlim          = extract<int>(ctrlist[13]);
    ctl.dt              = extract<double>(ctrlist[14]);
    // thermostat
    ctl.ntt             = extract<int>(ctrlist[15]);
    object obj_temp0    = extract<object>(ctrlist[16]);
    ctl.temp0           = (double *)(PyArray_DATA(obj_temp0.ptr()));
    object obj_tempi    = extract<object>(ctrlist[17]);
    ctl.tempi           = (double *)(PyArray_DATA(obj_tempi.ptr()));
    object obj_duration = extract<object>(ctrlist[18]);
    ctl.duration        = (int *)(PyArray_DATA(obj_duration.ptr()));
    ctl.ntemp           = extract<int>(ctrlist[19]);
    ctl.tautp           = extract<double>(ctrlist[20]);
    top.iCircular       = extract<int>(ctrlist[21]);
    ctl.verlet_cut      = extract<double>(ctrlist[22]);
    ctl.nVerlet         = extract<int>(ctrlist[23]);
    ctl.dpr             = extract<int>(ctrlist[24]);
    ctl.iCollide        = extract<int>(ctrlist[25]);
    printf("Done: extract control\n");

//---- Reading force field parameters ---------------------------
    //----------------------------------------------------
    printf("Start: extract force field parms\n");
    top.sigma           = extract<double>(prmlist[0]);
    top.epsilon         = extract<double>(prmlist[1]);
    top.vcut            = extract<double>(prmlist[2]);
    top.bphm_eq         = extract<double>(prmlist[3]);
    top.bphm_k          = extract<double>(prmlist[4]);
    top.bphm_p          = extract<double>(prmlist[5]);
    top.angle_eq        = extract<double>(prmlist[6]);
    top.angle_k         = extract<double>(prmlist[7]);
    top.angle_p         = extract<double>(prmlist[8]);
    top.bond_eq         = extract<double>(prmlist[9]);
    top.bond_k          = extract<double>(prmlist[10]);
    top.charge          = extract<double>(prmlist[11]);
    top.ecut            = extract<double>(prmlist[12]);
    top.len_Debye       = extract<double>(prmlist[13]);
    top.min_loop        = extract<int>(prmlist[14]);
    top.minbpFLAG       = extract<int>(prmlist[15]);
    top.minStembp       = extract<int>(prmlist[16]);
    top.bulge_p         = extract<double>(prmlist[17]);
    top.bcut            = extract<double>(prmlist[18]);
    top.stk_eq          = extract<double>(prmlist[19]);
    object obj_stk_k    = extract<object>(prmlist[20]);
    top.stk_k           = (double *)(PyArray_DATA(obj_stk_k.ptr()));
    top.stk_ks          = extract<double>(prmlist[21]);
    top.stk_p           = extract<double>(prmlist[22]);
    object obj_bpFLAG   = extract<object>(prmlist[23]);
    top.bpFLAG          = (int *)(PyArray_DATA(obj_bpFLAG.ptr()));
    object obj_nHB      = extract<object>(prmlist[24]);
    top.nHB             = (int *)(PyArray_DATA(obj_nHB.ptr()));
    top.f_scale         = extract<double>(prmlist[25]);

    //--------------------------------------------------------
    top.rcut            = ctl.rcut;
    top.stretch_force   = ctl.stretch_force;
    top.stretch_spring  = ctl.stretch_spring;
    top.restraint_wt   = ctl.restraint_wt;
    top.dist_init       = ctl.dist_init;
    top.dist_final      = ctl.dist_final;
    //top.ismd  = ctl.ismd;
    eng.W_smd = 0.0;
    

    //------------------------------------------------------------------------------------------------------
    // initialize restaint variables
    if (ctl.ntr or ctl.dpr){
        printf("Start: extract restraint\n");
        //------------- ntr ---------------------------------------------
        top.ntr_nBP      = extract<int>(ntrlist[0]);
        object obj_ntrBplist =    extract<object>(ntrlist[1]);
        top.ntrBplist  =    (int *)(PyArray_DATA(obj_ntrBplist.ptr()));
        //------------- dpr ---------------------------------------------
        top.dpr_nBP      = extract<int>(dprlist[0]);
        object obj_dprBplist =    extract<object>(dprlist[1]);
        top.dprBplist  =    (int *)(PyArray_DATA(obj_dprBplist.ptr()));
        printf("Done: extract restraint\n");
        n = 0;
        printf("ntr restriants:\n");
        for (i = 0; i < top.ntr_nBP; i++){
            printf("%d-%d  ", top.ntrBplist[n], top.ntrBplist[n+1]);
            n += 2;
        }
        printf("dpr restriants:\n");
        for (i = 0; i < top.dpr_nBP; i++){
            printf("%d-%d  ", top.dprBplist[n], top.dprBplist[n+1]);
            n += 2;
        }
        // multiply 3 to index
        ///for (i = 0; i < top.nrBp * 2; i++){
        ///    top.ntrBplist[i] *= 3;
        ///}
        /*
        for (i = 0; i < top.nchain; i++){
            n = i*top.maxHelix*2;
            printf("chain-%d has %d Basepairs: \n", i, top.nrBp[i]);
            for (j = 0; j < top.nrBp[i]; j++){
                // get atom index and force field parm index
                printf("%d-%d  ", top.ntrBplist[n], top.ntrBplist[n+1]);
                n += 2;
            }
            printf("\n");
        }*/
    }
    else{
        printf("Non-retrained MD!\n");
    }
    

    printf("the interface is okay!\n");
tim.t_Extract += punchClock(tim.tt);
    
    // output file !!!
    FILE *fp = NULL;
    FILE *fp_smd = NULL;
    fp = fopen(fnOut.c_str(), "w");

    // Time log
    time_t timep;        
    time (&timep); 
    printf("Current time: %s", asctime(gmtime(&timep)));
    fprintf(fp, "\n*******************************************************************************\n");
    fprintf(fp, "  Time log (GMT): %s", asctime(gmtime(&timep)));
    fprintf(fp, "*******************************************************************************\n");

    //---- steered MD log file ------------------------------------------------------------------
    if (ctl.ismd) {
        fp_smd = fopen("dist.RST", "w");
    }
    //---------------  C O N T R O L S ---------------------------------------------------------
    fprintf(fp, "\n*******************************************************************************\n");
    fprintf(fp, "  PlanarFold's CONTROL VARIABLES \n");
    fprintf(fp, "*******************************************************************************\n");
    fprintf(fp, "ntpr=%-d   ntwx=%-d  ntr=%-d\n", ctl.ntpr, ctl.ntwx, ctl.ntr);
    fprintf(fp, "ismd=%-d   dumpfreq=%-d\n", ctl.ismd, ctl.dumpfreq);
    fprintf(fp, "stretch_force=%-lg   stretch_spring=%-lg\n", ctl.stretch_force, ctl.stretch_spring);
    fprintf(fp, "dist_init=%-lg   dist_final=%-lg\n", ctl.dist_init, ctl.dist_final);
    fprintf(fp, "restraint_wt=%-lg   restraint_eq=%-lg\n", ctl.restraint_wt, ctl.restraint_eq);
    fprintf(fp, "rcut=%-lg   nscm=%-d\n", ctl.rcut, ctl.nscm);
    fprintf(fp, "nstlim=%-d   dt=%-lg\n", ctl.nstlim, ctl.dt);
    fprintf(fp, "ntt=%-d   ntemp=%-d   tautp=%-lg\n", ctl.ntt, ctl.ntemp, ctl.tautp);
    for (i = 0; i < ctl.ntemp; i++){
        fprintf(fp, "temp0[%d]=%-lg  tempi[%d]=%-lg  duration[%d]=%-d\n",i,ctl.temp0[i],i,ctl.tempi[i],i,ctl.duration[i]);
    }
    fprintf(fp, "iCircular=%-d\n", top.iCircular);
    fprintf(fp, "verlet_cut=%-lg   nVerlet=%-d\n", ctl.verlet_cut, ctl.nVerlet);
    fprintf(fp, "dpr=%-d\n", ctl.dpr);
    fprintf(fp, "iCollide=%-d\n", ctl.iCollide);


    //----------------- F O R C E   F I E L D ------------------------------------------------
    fprintf(fp, "\n*******************************************************************************\n");
    fprintf(fp, "  PlanarFold's FORCE FIELD PARAMETERS \n");
    fprintf(fp, "*******************************************************************************\n");
    fprintf(fp, "bond_eq=%-5lg     bond_k=%-5lg\n", top.bond_eq, top.bond_k);
    fprintf(fp, "sigma=%-5lg       epsilon=%-5lg   vcut=%-5lg\n", top.sigma, top.epsilon, top.vcut);
    fprintf(fp, "bphm_eq=%-5lg     bphm_k=%-5lg    bphm_p=%-5lg\n", top.bphm_eq, top.bphm_k, top.bphm_p);
    fprintf(fp, "min_loop=%-5d    minbpFLAG=%-5d minStembp=%-5d\n", top.min_loop, top.minbpFLAG, top.minStembp);
    fprintf(fp, "bulge_p=%-5lg     bcut=%-5lg\n", top.bulge_p, top.bcut);
    fprintf(fp, "angle_eq=%-5lg  angle_k=%-5lg   angle_p=%-5lg\n", top.angle_eq, top.angle_k, top.angle_p);
    fprintf(fp, "charge=%-5lg      ecut=%-5lg      len_Debye=%-5lg\n", top.charge, top.ecut, top.len_Debye);
    fprintf(fp, "stk_eq=%-5lg     stk_ks=%-5lg    stk_p=%-5lg\n", top.stk_eq, top.stk_ks, top.stk_p);
    fprintf(fp, "stk_k:\n");
    for (i = 0; i < 6; i++){
        for (j = 0; j < 6; j++){
            fprintf(fp, "%5.2lf ", top.stk_k[i*6+j]);
        }
        fprintf(fp, "\n");
    }
    for (i = 0; i < 2; i++){
        for (j = 0; j < 5; j++){
            fprintf(fp, "bpFLAG[%-d]=%-d  ", i*5+j, top.bpFLAG[i*5+j]);
        }
        fprintf(fp, "\n");
    }
    for (i = 0; i < 2; i++){
        for (j = 0; j < 5; j++){
            fprintf(fp, "nHB[%-d]=%-d  ", i*5+j, top.nHB[i*5+j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "f_scale=%-5lg\n", top.f_scale);

    //----------------- F O R C E   F I E L D ------------------------------------------------
    fprintf(fp, "\n*******************************************************************************\n");
    fprintf(fp, "  PlanarFold's RESTRAINTS / CONFINEMENTS \n");
    fprintf(fp, "*******************************************************************************\n");
    fprintf(fp, "ntr_nBP=%-d\n", top.ntr_nBP);
    if (ctl.ntr > 0){
        k = 1;
        for (j = 0; j < top.ntr_nBP; j++){
            fprintf(fp, "%3d-%-3d ", top.ntrBplist[j*2], top.ntrBplist[j*2+1]);
            if(k%10==0){
                fprintf(fp, "\n");
            }
            k++;
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "dpr_nBP=%-d\n", top.dpr_nBP);
    if (ctl.dpr > 0){
        k = 1;
        for (j = 0; j < top.dpr_nBP; j++){
            fprintf(fp, "%3d-%-3d ", top.dprBplist[j*2], top.dprBplist[j*2+1]);
            if(k%10==0){
                fprintf(fp, "\n");
            }
            k++;
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n*******************************************************************************\n");
    fprintf(fp, "  PlanarFold's ENERGIES \n");
    fprintf(fp, "*******************************************************************************\n");

//-------------------------------------------------------------------------
    // preparation work facilitating verlet integration loop
    prepare(&ctl, &top, &eng, &sdt, &tim, fp, 0);
	printf("the prepare() is done!\n");
    // cout<<"sizeof(sdt.acl[0]):"<<sizeof(sdt.acl[0])<<endl;

    // creat restart file (e.g. run001.rst)
    strcpy(ncf.Rst_filename, "");
    //strcat(ncf.Rst_filename, "run.rst");
    strcat(ncf.Rst_filename, fnRst.c_str());
    printf("ncf.Rst_filename: %s\n", ncf.Rst_filename);
    printf("Natom=%d\n", sdt.natom);
    createRst(sdt.natom, &ncf);
    printf("createRst is okay\n");

    // creat trajectory file (e.g. run001.nc)
    strcpy(ncf.Trj_filename, "");
    //strcat(ncf.Trj_filename, "run.nc");
    strcat(ncf.Trj_filename, fnTrj.c_str());
    printf("ncf.Trj_filename: %s\n", ncf.Trj_filename);
    nframe = (int) (ctl.nstlim / ctl.ntwx);     // e.g. nframe= 500,000/500 =1000
    createTrj(sdt.natom, nframe, &ncf);

    // loop for force and energy calculation
    int ntwx, ntpr;
    float *crdF;
    nstep = ctl.nstlim;
    ntwx = ctl.ntwx;
    ntpr = ctl.ntpr;
    // to store single precision coordinates for writeTrj
    crdF = (float*) malloc(3 * sdt.natom * sizeof(float));
    printf("before verlet is okay\n");

    // the thermostat: heating or cooling or equiling?
    int dura, duratot;
    double* dtemp;
    if (ctl.ntt==1){
        dtemp = (double*) malloc(ctl.ntemp * sizeof(double));
        duratot = 0;
        for(i = 0; i < ctl.ntemp; i++){
            dtemp[i] = (ctl.temp0[i] - ctl.tempi[i]) / ctl.duration[i];
            duratot += ctl.duration[i];
    }
    }
    // constant-velocity SMD
    double smd_vel;
    smd_vel =  (ctl.dist_final - ctl.dist_init) / nstep;

    // Preparation works are all done, now perform dynamics loop
tim.t_Prepare += punchClock(tim.tt);
    printf("verlet loop ... ... \n");
    for(i = 1; i < nstep+1; i++){
        /*
        printf("//-------------------------------------------------------------------\n");
        printf("//------------------------- L O O P %1d------------------------------\n", i);
        printf("//-------------------------------------------------------------------\n");
        */
        // current temperature
        if (ctl.ntt==1){
            dura = (i-1) % duratot;
            //------------------------------------------------------------------------------------------------------
            // to re-set the Ttarget to prevent "double precision shift" that will fail "restarting"!!!
            if (dura == 0){
                sdt.Ttarget = ctl.tempi[0];
            }
            //------------------------------------------------------------------------------------------------------
            for (j = 0; j < ctl.ntemp; j++){
                if (dura < ctl.duration[j]){
                    sdt.Ttarget += dtemp[j];
                    break;
                }
                else{
                    dura -= ctl.duration[j];
                }
            }
        }
        // update the equil value for constant-velocity SMD
        if (ctl.ismd == 2) {
            if (i % ctl.dumpfreq == 0){
                sdt.spring_len = ctl.dist_init + i * smd_vel;
            }
        }
        verlet(&ctl, &top, &eng, &sdt, &ncf, &tim, i, fp);
        // Write trajectory file
        if (i % ntwx == 0){  //i>0 cause we don't wish crd[0] to be stored
            writeTrj((int)(i/ntwx)-1, &sdt, &ncf, crdF);
tim.t_Write += punchClock(tim.tt);
        }
        // Write SMD-related outputs
        if (ctl.ismd) {
            if (i % ctl.dumpfreq == 0){
                fprintf(fp_smd, "%10.5lf %10.5lf %10.5lf %10.5lf\n", sdt.spring_len, sdt.dist_end, top.stretch_force, eng.W_smd);
            }
        }
        // Write energies
        if (i % ntpr == 0){  //i>0 cause we don't wish crd[0] to be stored
            //printf("//------------------------- L O O P %1d-----------------------------------\n", i);
            /*
            printf("Bond:        %8.4lf    Vdw:              %8.4lf    EPtot:      %8.4lf\n", eng.E_bond, eng.E_vdw, eng.Epot);
            printf("BasePair:    %8.4lf    Angle:            %8.4lf    Stacking:   %8.4lf\n", eng.E_basepair, eng.E_angl, eng.E_stacking);
            printf("EKtot:   %12.4lf    Temperature:  %12.4lf    Etot:   %12.4lf\n", eng.Ekin, sdt.T, eng.Etot);
            */
            //outputfiles
            fprintf(fp, "//-----------------------------------------------------------------------------\n");
            fprintf(fp, "Frame:   %12d   Temperature:    %12.4lf   Restraint:  %12.4lf\n", i, sdt.T, eng.E_restraint);
            fprintf(fp, "Bond:    %12.4lf   Angle:          %12.4lf   Stacking:   %12.4lf\n", eng.E_bond, eng.E_angl, eng.E_stacking);
            fprintf(fp, "Vdw:     %12.4lf   Elec:           %12.4lf   BasePair:   %12.4lf\n", eng.E_vdw, eng.E_elec, eng.E_basepair);
            fprintf(fp, "Ektot:   %12.4lf   Eptot:          %12.4lf   Etot:       %12.4lf\n", eng.Ekin, eng.Epot, eng.Etot);
            //fprintf(fp, "Ttarget: %12.4lf\n", sdt.Ttarget);
            fprintf(fp, "!!BasePairs: %d\n", top.nBP);
            k = 1;
            for (j = 0; j < top.nBP; j++){
                fprintf(fp, "%3d-%-3d ", top.aBPlist5[j], top.aBPlist3[j]);
                if(k%10==0){
                    fprintf(fp, "\n");
                }
                k++;
            }
            if((k-1)%10!=0){
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n");
            // test
            if (ctl.ntr == 2){
                fprintf(fp, "!!targeted BP: %d\n", top.n_tgtBP);
                k = 1;
                for (j = 0; j < top.n_tgtBP; j++){
                    fprintf(fp, "%3d-%-3d ", top.tgtBplist[j*2], top.tgtBplist[j*2+1]);
                    if(k%10==0){
                        fprintf(fp, "\n");
                    }
                    k++;
                }
                fprintf(fp, "\n");
            }
            // test -----------------------------------------------------------------------
            /*
            fprintf(fp, "sdt.verlet_natom and sdt.verlet_list\n");
            for (j = 0; j < sdt.natom-2; j++){
                fprintf(fp, "%d(%d):", j, sdt.verlet_natom[j]);
                for (k = 0; k < sdt.verlet_natom[j]; k++){
                    fprintf(fp, "%d \t", sdt.verlet_list[j][k]);
                }
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n");
            */

            // test
            /*
            fprintf(fp, "coordinates and volocities\n");
            for (j = 0; j < sdt.natom; j++){
                fprintf(fp, "%d, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf\n", j, sdt.crd[j*3], sdt.crd[j*3+1], sdt.vel[j*3], sdt.vel[j*3+1], sdt.acl[j*3], sdt.acl[j*3+1]);
            }
            fprintf(fp, "\n");
            */
            //
            /*
            k = 1;
            fprintf(fp, "!!Stackings: %d\n", top.nSTk);
            for (j = 0; j < top.nSTk; j++){
                fprintf(fp, "%d ", top.stklist[j]);
                if(k%10==0){
                    fprintf(fp, "\n");
                }
                k++;
            }
            fprintf(fp, "\n");
            */
tim.t_Write += punchClock(tim.tt);
        }
    }
    printf("the verlet loop is done!\n");

    // 1.write restart file, and 2.close rst/nc file
    if (i == nstep+1){
        //writeTrj(nframe-1, &sdt, &ncf, crdF);
        writeRst(nframe, &sdt, &ncf);
        tim.t_Write += punchClock(tim.tt);
    }
    free(crdF);
    printf("Trj and Nc files are written!\n");
    closeRstTrj(&ncf);
    // deallocate the pointers allocated in 'prepare()'
	printf("the freeC() is done!\n");

    // stop the time logging
    StopTimer(&tim, fp);
    fclose(fp);
    if (ctl.ntr == 3) {
        fclose(fp_smd);
    }
    printf("output file is written and closed!\n");
}

BOOST_PYTHON_MODULE(dsflyext)
{
    def("run_md", run_md);
}


