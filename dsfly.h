#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sys/time.h>
#define PI ((double)3.1415926535898)
// The universal gas constant in kcal/mol-K
#define GASCNST ((double)1.987216e-3)
// The maximum excluded-atom-list length
#define MAXNEX ((int)25)


//-----------------------------------S T R U C T-------------------------------
struct OutputNetCDFfiles {
    // Trajectory file
    int ncid_Trj;
    char Trj_filename[10];    // "runNNN.nc", N=0~9, (Trj_filename[10] is for '\0')
    int varid_time_Trj;      // The variable ID for 'time' of trajectory file
    int varid_spatial_Trj;   // The variable ID for 'spatial' of trajectory file
    int varid_crd_Trj;       // The variable ID for 'coordinates' of trajectory file
    int varid_cellSpatial_Trj;   // The variable ID for 'cell_spatial' of trajectory file
    int varid_cellAngular_Trj;   // The variable ID for 'cell_angular' of trajectory file
    int varid_cellLengths_Trj;   // The variable ID for 'cell_lengths' of trajectory file
    int varid_cellAngles_Trj;    // The variable ID for 'cell_angles' of trajectory file

    // Restart file
    int ncid_Rst;
    char Rst_filename[11];    // "runNNN.rst", N=0~9, (Crd_filename[11] is for '\0')
    int varid_time_Rst;     // The variable ID for 'time' of restart file
    int varid_spatial_Rst;  // The variable ID for 'spatial' of restart file
    int varid_crd_Rst;      // The variable ID for 'coordinates' of restart file
    int varid_vel_Rst;      // The variable ID for 'velocities' of restart file
    int varid_cellSpatial_Rst;  // The variable ID for 'cell spatial' of restart file
    int varid_cellAngular_Rst;  // The variable ID for 'cell_angular' of restart file
    int varid_cellLengths_Rst;  // The variable ID for 'cell_lengths' of restart file
    int varid_cellAngles_Rst;   // The variable ID for 'cell_angles' of restart file
};
typedef struct OutputNetCDFfiles ncfile;

struct States {
    int natom;          // The number of atoms in the system
    int ndf;            // The number of degrees of freedom in the system
    double Ttarget;     // Target bath temperature
    double T;           // current instantaneous temperature
    double* fcrd;       // The fractional coordinates of the atoms
    double* crd;        // The (current) coordinates of the atoms
    double* vel;        // The (current) velocities of the atoms
    double* acl;        // The (current) acceleration on atoms
    double box[6];      // The simulation box lengths and angles
    int** NotSearched;        // The Flag denoting whether this nucleotide is Base-Pair-saturate
    double** cMap;      // Contact map
    double** triMap;    // triangularized Contact map
    double** sumMap;    // triangularized Contact map
    // triplex exclusion
    int** col_iblock;   // the index of the bolck (column) of this cmap[i,j]
    int** row_iblock;   // the index of the bolck (row) of this cmap[i,j]
    int* max_iblock;    // the index of the maximum of the bolck (row or column, symmetric)
    //
    int* stlist;        // stem Start-points list in Contact map
    int* ndlist;        // stem End-points list in Contact map
    // SMD-realted
    double smd_frc;     // the force of the stretching force
    double dist_end;    // the distance between chain termini (where force has been applied to)
    double spring_len;  // the equilibrium value of spring
    // verlet-list related
    int** verlet_list;
    int* verlet_natom;
    long int rndcon;    // The random number counter
};
typedef struct States state;

struct Topology{
    int* pointers;      // POINTERS
    double* charges;     // CHARGE
    double* mass;       // MASS
    int* iac;           // ATOM_TYPE_INDEX
    int* ico;           // NONBONDED_PARM_INDEX
    int* inb;           // EXCLUDED_ATOMS_LIST
    int* nblist;        // The residues need calculation of vdW and Elec
    int nNBres;         // number of excluded residues
    int nbres;          // number of remaining residues
    // Force field parameters read from the input file
    double sigma;       // vdW interaction
    double epsilon;     // vdW interaction
    double vcut;        // vdW cutoff
    double bphm_eq;     // BP-harmonic equilibrium value
    double bphm_k;      // BP-harmonic force constant
    double bphm_p;      // BP-harmonic score bonus
    double angle_eq;    // Angle equilibrium value
    double angle_k;     // Angle force constant
    double angle_p;     // Angle score bonus
    double bond_eq;     // Bond equilibrium value
    double bond_k;      // Bond force constant
    double charge;      // Net charge
    double ecut;        // Electrostatic cutoff
    double len_Debye;   // Debye length in electrostatic 
    int min_loop;       // minimal loop size in dynamic programming. (!! Maximum: 5)
    int minbpFLAG;      // minimum of bpFLAG to calculate the cMap
    int minStembp;      // the minimum bp allowed in a 'stem'
    double bulge_p;     // penalty of bulge in dynamic programming
    double bcut;        // BP-harmonic cutoff
    double stk_eq;      // Stacking equilibrium value
    double* stk_k;      // Stacking force constant (Turner's 36 parms)
    double stk_ks;      // Stacking force constant scaling factor
    double stk_p;       // Stacking score bonus
    int* bpFLAG;        // The number of HB for each
    int* nHB;
    double f_scale;     // Uniform scaling factor for all potentials
    // 
    int* anglist;       // ANGLES_WITHOUT_HYDROGEN
    int nangle;
    double* inv_mass;   // The inverse atomic masses                                    // need to be supplemented before gradient()
    double* cn1;        // LENNARD_JONES_ACOEF
    double* fcn1;       // LENNARD_JONES_ACOEF used for force calculation(*12)
    double* cn2;        // LENNARD_JONES_BCOEF
    double* fcn2;       // LENNARD_JONES_BCOEF used for force calculation(*6)
    double restraint_wt;    // Force constant for ntr=1 use
    double restraint_eq;    // Equil value for ntr=1 use
    double ljA;
    double ljB;
    double fljA;
    double fljB;
    //
    int* BPlist5;
    int* BPlist3;
    int* aBPlist5;
    int* aBPlist3;
    //
    int** BP_partners;
    int* nBP_partner;
    //
    double* cmap_max_row;
    int nBP;
    int nStem;
    int* Stem;
    int* aStem;
    int* isTermi;
    int* BPlistTp;
    int* stklist;
    int nSTk;
    int* isStem;
    // restrainted MD, ntr
    int* ntrBplist;
    int ntr_nBP;
    double rcut;
    // restrainted MD, dpr
    int* dprBplist;
    int dpr_nBP;
    // targeted MD
    int* tgtBplist;
    int* notTargeted;
    int n_tgtBP;
    int* orient;
    //SMD-related
    double stretch_force;
    double dist_init;
    double dist_final;
    double MassTot;
    // circular RNA
    int iCircular;
    // steered MD
    //int ismd;
    double stretch_spring;  // the spring constant for c.v.SMD
};
typedef struct Topology prmtop;

struct Energy{
    double E_bond;          // Energy due to bond-stretching terms
    double E_angl;          // Energy due to angle-bending terms
    double E_dihe;          // Energy due to dihedral terms
    double E_vdw;           // Energy due to van der Waarl terms
    double E_elec;          // Energy due to electrostatic terms
    double E_stacking;      // Energy due to base-stacking terms
    double E_basepair;      // Energy due to base-pairing terms
    double E_BPharmo;       // Energy due to Bp-harmonic terms
    double E_restraint;     // Energy due to restraint terms
    double Epot;            // Total potential energy
    double Ekin;            // Total kinetic energy
    double Etot;            // Total energy of the system
    // SMD-realted
    double E_smd;
    double W_smd;        // the work summed along the trajectory
};
typedef struct Energy energy; 

struct Mdcontrol{
    int ntpr;       // Writing frequency of output diagnostics file (mdout)
    int ntwx;       // Writing frequency of coordinate trajectory file (mdcrd)
    //
    int ntr;        // Flag for restraining specified atoms in Cartesian
    int dpr;        // Flag for confined-MD: dpr=1, no DP at all; dp=2, only neighboring bps in DP
    // SMD related
    int ismd;       // option for Steered MD
    int dumpfreq;
    double stretch_force;   // the stretching force for c.f.SMD, unit in kcal/mol/A
    double stretch_spring;  // the spring constant for c.v.SMD
    double dist_init;
    double dist_final;
    double restraint_wt;    // Force constant for ntr=1 use
    double restraint_eq;    // Equil value for ntr=1 use
    double rcut;
    // 
    int nscm;      // Frequency of removing translation and rotation of COM
    int nstlim;     // Number of MD-steps to be performed
    double dt;      // Time steps (ps)
    // temperature control
    int ntt;        // Thermostat type, 0/1 for no/Berendsen Thermostat
    double* temp0;   // Reference temperature for ntt>0
    double* tempi;   // Initial temperature
    int* duration;   // The duration of thermostat
    int ntemp;     // the section of different thermostating stage
    double tautp;   // Time constant for temperature fluctuations (in ps)
    // verlet list updating control
    double verlet_cut;     // the verlet list cutoff (~ 30A would be fine)
    int nVerlet;    // the step to update verlet list (1000 would be fine)
    // the random collision with water to avoid sine wave (for very long non-pairing RNA)
    int iCollide;   // to asign random collision of water with residues
};
typedef struct Mdcontrol control;

struct TimeLog{
    // time Unit : second/s !!!
    struct timeval tt[2];   // time of start point tt[0], and end point tt[1] (merge mdgx's tti and ttf together for easier variable)
    struct timeval t_Start;  // time point data passed from python into C/C++.
    struct timeval t_End;    // time point when all C/C++ data are freed and missons are finished.
    double t_Extract;       // time of the extraction of input informations from python to C/C++ using boost library
    double t_Prepare;       // time of all preparation works
    double t_Verlet;        // time of verlet integration loop updating vel and crd
    double t_Bond;          // time of calculating Bond interaction
    double t_Angle;
    double t_Stacking;
    double t_VDW;           // time of calculating van der Waals interaction
    double t_ELEC;          // time of calculating electrostatic interaction
    double t_BasePair;      // time of calculating base-pairing interaction
    double t_Restraint;     // time of calculating restraint
    double t_Thermostat;    // time of thermostating
    double t_Barostat;      // time of barostating
    double t_Write;         // time of writing output files, including rst, nc, and mdout.
    double t_BP_DP;   
    double t_BP_others;
    double t_BP_HB;
    double t_BP_exTri;
    double t_BP_exPK;
    double t_BP_Angs;
    //
    int* Memory_check;
};
typedef struct TimeLog timer;

//------------------------------F U N C T I O N--------------------------------
void InitAll(control *ctl, prmtop *top, energy *eng, state *sdt, ncfile *ncf);
void InitTimer(timer *tim);
void StopTimer(timer *tim, FILE *fp);
double punchClock(struct timeval tt[]);

double perfc(double x);
double PerfcFast(double x);
void CrossP(double* p, double* q, double* cr);
void UnitVector3(double* V);
double Dot3(double* V1, double* V2);

void prepare(control *ctl, prmtop *top, energy *eng,
             state *sdt, timer *tim, FILE *fp, int nstep);
void verlet(control *ctl, prmtop *top, energy *eng,
        state *sdt, ncfile *ncf, timer *tim, int nstep, FILE *fp);
void gradient(energy *eng, prmtop *top, state *sdt,
            timer *tim, control *ctl, FILE *fp, int nstep);

void createRst(int natom, ncfile *ncf);
void createTrj(int natom, int nframe, ncfile *ncf);
void writeTrj(int currFrame, state *sdt, ncfile *ncf, float *crdF);
void writeRst(int nframe, state *sdt, ncfile *ncf);
void closeRstTrj(ncfile *ncf);
