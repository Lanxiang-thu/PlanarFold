//#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include "dsfly.h"
// This is referencing mdgx's Timingd.c

//-----------------------------------------------------------------------------
// punchClock: record the time just after finishing some piece of code.
//             Before exiting, this function automatically sets the initial
//             time to the final time, so that multiple calls to punchClock()
//             can be used to record a number of successive time intervals.
//                                                                      
// Arguments:                                                           
//   tt:    the timer structs stored the start and end point of a time interval 
//-----------------------------------------------------------------------------
double punchClock(struct timeval tt[]){
    double dt;

    gettimeofday(&tt[1], NULL);
    dt = tt[1].tv_sec - tt[0].tv_sec +
        (1.0e-6)*(tt[1].tv_usec - tt[0].tv_usec);
    tt[0] = tt[1];

    return dt;
}

//-----------------------------------------------------------------------------
// InitTimer: initialize a structure to keep timing information                                         
//-----------------------------------------------------------------------------
void InitTimer(timer *tim){
    tim->t_Extract = 0.0;       // time of the extraction of input informations from python to C/C++ using boost library
    tim->t_Prepare = 0.0;       // time of all preparation works
    tim->t_Verlet = 0.0;        // time of verlet integration loop updating vel and crd. (just the integration time)
    tim->t_Bond = 0.0;          // time of calculating Bond interaction
    tim->t_Angle = 0.0;
    tim->t_Stacking = 0.0; 
    tim->t_VDW = 0.0;           // time of calculating van der Waals interaction
    tim->t_ELEC = 0.0;          // time of calculating electrostatic interaction
    tim->t_BasePair = 0.0;      // time of calculating base-pairing interaction
    tim->t_BP_DP = 0.0;  
    tim->t_BP_exTri = 0.0;
    tim->t_BP_exPK = 0.0;
    tim->t_BP_Angs = 0.0;
    tim->t_BP_others = 0.0;
    tim->t_BP_cmap = 0.0;
    tim->t_BP_blocks = 0.0;
    tim->t_Thermostat = 0.0;    // time of thermostating. (including the EKtot computation)
    tim->t_Restraint = 0.0;     // time of calculating restraint
    tim->t_Barostat = 0.0;      // time of barostating
    tim->t_Write = 0.0;         // time of writing output files, including rst, nc, and mdout.
    // record the time just before beginning any piece of C/C++ code.
    gettimeofday(&tim->t_Start, NULL);
    tim->tt[0] = tim->t_Start;
}

//-----------------------------------------------------------------------------
// StopTimer:                                      
//-----------------------------------------------------------------------------
void StopTimer(timer *tim, FILE *fp){
    double dtot;

    gettimeofday(&tim->t_End, NULL);
    dtot = tim->t_End.tv_sec - tim->t_Start.tv_sec +
        (1.0e-6)*(tim->t_End.tv_usec - tim->t_Start.tv_usec);

    tim->t_BasePair = tim->t_BP_DP + tim->t_BP_others + tim->t_BP_HB
                    + tim->t_BP_exTri + tim->t_BP_exPK + tim->t_BP_Angs
		    + tim->t_BP_cmap + tim->t_BP_blocks;

    printf("This is a memery check");
    tim->Memory_check = (int *)malloc(10*sizeof(int));
    printf("No memery corruption !");
    
    // exhibit time
    printf("\n DRfold's TIMINGS \n");
    printf(" Segment       |     Time(s)  /  Percentage\n");
    printf(" --------------+---------------------------\n");
    printf(" Extract       | %10.2lf       %6.2lf\n", tim->t_Extract, 100 * tim->t_Extract/dtot);
    printf(" Prepare       | %10.2lf       %6.2lf\n", tim->t_Prepare, 100 * tim->t_Prepare/dtot);
    printf(" Verlet        | %10.2lf       %6.2lf\n", tim->t_Verlet, 100 * tim->t_Verlet/dtot);
    printf(" Bond          | %10.2lf       %6.2lf\n", tim->t_Bond, 100 * tim->t_Bond/dtot);
    printf(" Angle         | %10.2lf       %6.2lf\n", tim->t_Angle, 100 * tim->t_Angle/dtot);
    printf(" Stacking      | %10.2lf       %6.2lf\n", tim->t_Stacking, 100 * tim->t_Stacking/dtot);
    printf(" t_VDW         | %10.2lf       %6.2lf\n", tim->t_VDW, 100 * tim->t_VDW/dtot);
    printf(" t_ELEC        | %10.2lf       %6.2lf\n", tim->t_ELEC, 100 * tim->t_ELEC/dtot);
    printf(" t_BasePair    | %10.2lf       %6.2lf\n", tim->t_BasePair, 100 * tim->t_BasePair/dtot);
    printf("    t_BP_DP       | %10.2lf       %6.2lf\n", tim->t_BP_DP, 100 * tim->t_BP_DP/dtot);
    printf("    t_BP_others   | %10.2lf       %6.2lf\n", tim->t_BP_others, 100 * tim->t_BP_others/dtot);
    printf("    t_BP_cmap     | %10.2lf       %6.2lf\n", tim->t_BP_cmap, 100 * tim->t_BP_cmap/dtot);
    printf("    t_BP_blocks   | %10.2lf       %6.2lf\n", tim->t_BP_blocks, 100 * tim->t_BP_blocks/dtot);
    printf("    t_BP_exTri    | %10.2lf       %6.2lf\n", tim->t_BP_exTri, 100 * tim->t_BP_exTri/dtot);
    printf("    t_BP_exPK     | %10.2lf       %6.2lf\n", tim->t_BP_exPK, 100 * tim->t_BP_exPK/dtot);
    printf("    t_BP_Angs     | %10.2lf       %6.2lf\n", tim->t_BP_Angs, 100 * tim->t_BP_Angs/dtot);
    printf(" t_Restraint   | %10.2lf       %6.2lf\n", tim->t_Restraint, 100 * tim->t_Restraint/dtot);
    printf(" Thermostat    | %10.2lf       %6.2lf\n", tim->t_Thermostat, 100 * tim->t_Thermostat/dtot);
    printf(" Write         | %10.2lf       %6.2lf\n", tim->t_Write, 100 * tim->t_Write/dtot);
    printf(" --------------+---------------------------\n");
    printf(" Total CPU Time: %10.2lf       %6.2lf\n", dtot, 100.0);
    
    // store time message
    fprintf(fp, "\n***********************************************************************\n");
    fprintf(fp, "  DRfold's TIMINGS \n");
    fprintf(fp, "***********************************************************************\n");
    fprintf(fp, " Segment       |     Time(s)  /  Percentage\n");
    fprintf(fp, " --------------+---------------------------\n");
    fprintf(fp, " Extract       | %10.2lf       %6.2lf\n", tim->t_Extract, 100 * tim->t_Extract/dtot);
    fprintf(fp, " Prepare       | %10.2lf       %6.2lf\n", tim->t_Prepare, 100 * tim->t_Prepare/dtot);
    fprintf(fp, " Verlet        | %10.2lf       %6.2lf\n", tim->t_Verlet, 100 * tim->t_Verlet/dtot);
    fprintf(fp, " Bond          | %10.2lf       %6.2lf\n", tim->t_Bond, 100 * tim->t_Bond/dtot);
    fprintf(fp, " Angle         | %10.2lf       %6.2lf\n", tim->t_Angle, 100 * tim->t_Angle/dtot);
    fprintf(fp, " Stacking      | %10.2lf       %6.2lf\n", tim->t_Stacking, 100 * tim->t_Stacking/dtot);
    fprintf(fp, " t_VDW         | %10.2lf       %6.2lf\n", tim->t_VDW, 100 * tim->t_VDW/dtot);
    fprintf(fp, " t_ELEC        | %10.2lf       %6.2lf\n", tim->t_ELEC, 100 * tim->t_ELEC/dtot);
    fprintf(fp, " t_BasePair    | %10.2lf       %6.2lf\n", tim->t_BasePair, 100 * tim->t_BasePair/dtot);
    fprintf(fp, "    t_BP_DP       | %10.2lf       %6.2lf\n", tim->t_BP_DP, 100 * tim->t_BP_DP/dtot);
    fprintf(fp, "    t_BP_others   | %10.2lf       %6.2lf\n", tim->t_BP_others, 100 * tim->t_BP_others/dtot);
    fprintf(fp, "    t_BP_exTri    | %10.2lf       %6.2lf\n", tim->t_BP_exTri, 100 * tim->t_BP_exTri/dtot);
    fprintf(fp, "    t_BP_exPK     | %10.2lf       %6.2lf\n", tim->t_BP_exPK, 100 * tim->t_BP_exPK/dtot);
    fprintf(fp, "    t_BP_Angs     | %10.2lf       %6.2lf\n", tim->t_BP_Angs, 100 * tim->t_BP_Angs/dtot);
    fprintf(fp, " t_Restraint   | %10.2lf       %6.2lf\n", tim->t_Restraint, 100 * tim->t_Restraint/dtot);
    fprintf(fp, " Thermostat    | %10.2lf       %6.2lf\n", tim->t_Thermostat, 100 * tim->t_Thermostat/dtot);
    fprintf(fp, " Write         | %10.2lf       %6.2lf\n", tim->t_Write, 100 * tim->t_Write/dtot);
    fprintf(fp, " --------------+---------------------------\n");
    fprintf(fp, " Total CPU Time: %10.2lf       %6.2lf\n", dtot, 100.0);
}
