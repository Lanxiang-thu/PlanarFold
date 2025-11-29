#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dsfly.h"
#define PI ((double)3.1415926535898)
#define MaxNest ((int)2)

//-----------------------------------------------------------------------------
// ebond: computes the force on two beads due to bond stretching.
// 
//-----------------------------------------------------------------------------
void ebond(prmtop *top, state *sdt, double* E_bond){
    int i, nBond;
    int ia, ib, n;
    double eq, k;
    double dx, dy, de;
    double r2, r, dl;
    double *eng = NULL; // bond energy
    double *deriv = NULL; // first derivative
    int iCircular = top->iCircular;

    //
    eng = E_bond;
    deriv = sdt->acl;

    // The number of bonds for a polymer chain: Natom-1
    nBond = sdt->natom - 1;
    n = 0;
    // All bonds adopt uniform parameters
    eq = top->bond_eq;
    k = top->bond_k;
    //printf("nBond=%d\n", nBond);
    for (i = 0; i < nBond; i++){
        // get atom index
        ia = n;
        ib = n + 3;
        n += 3;
        // Compute displacement
        dx = sdt->crd[ia] - sdt->crd[ib];
        dy = sdt->crd[ia+1] - sdt->crd[ib+1];
        r2 = dx*dx + dy*dy;
        r = sqrt(r2);
        dl = r - eq;
        //printf("ia=%d, ib=%d, r=%lf\n", ia, ib, r);
        // Accumulate the bond energy and first derivative
        *eng += k*dl*dl;
        de = 2.0*k*dl/r;
        //
        deriv[ia] += de*dx;
        deriv[ia+1] += de*dy;
        deriv[ib] -= de*dx;
        deriv[ib+1] -= de*dy;
    }
    /*
    // for the circular RNA: add a bond to res[0]-res[n-1]
    if (iCircular){
        ia = 0;
        ib = nBond*3;
        // Compute displacement
        dx = sdt->crd[ia] - sdt->crd[ib];
        dy = sdt->crd[ia+1] - sdt->crd[ib+1];
        r2 = dx*dx + dy*dy;
        r = sqrt(r2);
        dl = r - eq;
        // Accumulate the bond energy and first derivative
        *eng += k*dl*dl;
        de = 2.0*k*dl/r;
        //
        deriv[ia] += de*dx;
        deriv[ia+1] += de*dy;
        deriv[ib] -= de*dx;
        deriv[ib+1] -= de*dy;
    }
    */
}


//-----------------------------------------------------------------------------------------
// eNonBonded : computes all the non-bonded  interaction 
//              including vdW, electrostatic, and the preparation of contact map for BP
//-----------------------------------------------------------------------------------------
void eNonBonded(prmtop *top, state *sdt, double* E_vdw, double* E_elec, int dpr, FILE *fp, int nstep){
    int i, j, k;
    int natom, n, n2;
    int ia, ib, ia3, ib3;
    double r2, r, invr2, invr4, invr6, invr;    
    // cutoffs
    double rc2_v = top->vcut * top->vcut;   // van der Waals
    double rc2_e = top->ecut * top->ecut;   // electrostatics
    double rc2_b = top->bcut * top->bcut;   // Basepair
    
    // verlet list
    int **verlet_list = sdt->verlet_list;
    int *verlet_natom = sdt->verlet_natom;

    // force and acceleration
    double dx, dy, deE, deV, de;
    double *deriv = NULL; // first derivative
    double *engVdw = E_vdw, *engElec = E_elec;
    deriv = sdt->acl;
    natom = sdt->natom;

    // van der Waals
    double ljA, ljB, fljA, fljB;
    engVdw = E_vdw;
    ljA = top->ljA;
    ljB = top->ljB;
    fljA = top->fljA;
    fljB = top->fljB;

    // electrostatics
    double charge = top->charge;
    double qA, qB, qq;
    double L_Debye = top->len_Debye; //Debye_length in Angstrom
    engElec = E_elec;
    qA = charge;
    qB = charge;
    qq = qA*qB;

    // Basepair
    double **cMap= NULL;
    int iparm, idxa, idxb;
    int ntypes = top->pointers[1];
    int minLoop = top->min_loop;
    int **NotSearched = NULL;
    //  if dpr=1, cmap and DP are not needed at all.
    if (dpr-1){
        cMap = sdt->cMap;
        NotSearched = sdt->NotSearched;
        for (i = 0; i < natom; i++){
            for (j = 0; j < natom; j++){
                NotSearched[i][j] = 1;
                //NotSearched[j][i] = 1;
                cMap[i][j] = 0.0;       // necessary!
                //cMap[j][i] = 0.0;       // necessary!
            }
        }
    }
    //fprintf(fp, "check eNonBonded part-1\n");
    
    // ===== start to serach pair-wise interactions ===============================
    for (i = 0; i < natom-2; i++){
        ia3 = i*3;
        idxa = top->iac[i]-1;
        //fprintf(fp, "check eNonBonded part-2\n");
        for (j = 0; j < verlet_natom[i]; j++){
            ib = verlet_list[i][j];
            ib3 = ib * 3;
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            // Usually:  ecut(20.0) > bcut(14.0) > vcut(13.5)
            //----- 1. electrostatics ---------------------------------------------------------------
            if (r2 < rc2_e){
                r = sqrt(r2);
                //fprintf(fp, "ia=%d, ib=%d, r=%lf\t", i, ib, r);
                invr = 1.0 / r;
                deE = qq * invr * exp(-r/L_Debye);
                *engElec += deE;
                deE *= -(1.0/L_Debye + invr) * invr;
                //----- 2. Basepair ---------------------------------------------------------------
                // if dpr=0, cmap are calculated here;
                // if dpr=2, cmap calculated in eBasePair()
                // if dpr=1, DP and cmap are skipped
                if (dpr == 0){
                    if (r2 < rc2_b){
                        // ------------- 0 --- 1 --- 2 --- 3 --- 4 --- 5 --- 6 -----
                        // vdw:          i     -    [o     o     o     o     o
                        // BP(min=3):    i     -     x     x    [o     o     o
                        // BP(min=4):    i     -     x     x     x    [o     o
                        // ----------------------------------------------------------
                        // prohibit the residues within smallest apical loops
                        if ((ib-i) > minLoop){
                            idxb = top->iac[ib]-1;
                            iparm = top->ico[ntypes*idxa + idxb] -1;
                            cMap[i][ib] = top->bpFLAG[iparm] / r2;   // considering pairing type
                            cMap[ib][i] = cMap[i][ib];
                        }
                    }
                }
                //----- 3. van der Waals --------------------------------------------------------
                deV = 0.0;
                if (r2 < rc2_v){
                    invr2 = invr*invr;
                    invr4 = invr2*invr2;
                    invr6 = invr4*invr2;
                    *engVdw += invr6 * (ljA*invr6 - ljB);
                    deV = invr4 * invr4 * (fljB - fljA*invr6);
                }
                de = deE + deV;
                //fprintf(fp, "deE=%lf， deV=%lf\n", deE, deV);
                // accumulate the first derivative
                deriv[ia3] += de*dx;
                deriv[ia3+1] += de*dy;
                deriv[ib3] -= de*dx;
                deriv[ib3+1] -= de*dy;
            }
        }
    }
}


//-----------------------------------------------------------------------------
// evdw : computes the van der Waals pairwise interaction using Lennard-Jones
//       potential.
//-----------------------------------------------------------------------------
void evdw(prmtop *top, state *sdt, double* E_vdw){
    int i, j, k;
    int natom, n, n2;
    int ia, ib, ia3, ib3;
    double ljA, ljB, fljA, fljB;
    double dx, dy, dz, de;
    double r2, r, invr2, invr4, invr6, invr;
    double *eng = NULL;
    double *deriv = NULL; // first derivative
    double rc2 = top->vcut * top->vcut;
    //int iCircular = top->iCircular;
    int **verlet_list = sdt->verlet_list;
    int *verlet_natom = sdt->verlet_natom;

    //
    eng = E_vdw;
    deriv = sdt->acl;
    natom = sdt->natom;
    ljA = top->ljA;
    ljB = top->ljB;
    fljA = top->fljA;
    fljB = top->fljB;
    
    for (i = 0; i < natom-2; i++){
        ia3 = i*3;
        for (j = 0; j < verlet_natom[i]; j++){
            ib3 = verlet_list[i][j] * 3;
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > rc2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            r = sqrt(r2);
            invr = 1.0 / r;
            invr2 = invr*invr;
            invr4 = invr2*invr2;
            invr6 = invr4*invr2;
            // accumulate the bond energy and first derivative
            *eng += invr6 * (ljA*invr6 - ljB);
            de = invr4 * invr4 * (fljB - fljA*invr6);
            deriv[ia3] += de*dx;
            deriv[ia3+1] += de*dy;
            deriv[ib3] -= de*dx;
            deriv[ib3+1] -= de*dy;
        }
    }
    /*
    // for Excluded list 
    for (i = 0; i < top->nbres-1; i++){
        ia3 = top->nblist[i]*3;
        // j=i+2 as we exclude two adjacent nucleotides
        // G--A--U--C--'G'--A--C--G--U
        // o--o--o--x--'*'--x--o--o--o
        for (j = i+1; j < top->nbres; j++){
            ib3 = top->nblist[j]*3;
            if (ia3-ib3==3 || ib3-ia3==3){
                continue;
            }
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > rc2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            r = sqrt(r2);
            invr = 1.0 / r;
            invr2 = invr*invr;
            invr4 = invr2*invr2;
            invr6 = invr4*invr2;
            // accumulate the bond energy and first derivative
            *eng += invr6 * (ljA*invr6 - ljB);
            de = invr4 * invr4 * (fljB - fljA*invr6);
            deriv[ia3] += de*dx;
            deriv[ia3+1] += de*dy;
            deriv[ib3] -= de*dx;
            deriv[ib3+1] -= de*dy;
        }
    }
    */
    /*
    // for the circular RNA (not compatible with excluded list yet)
    for (i = 0; i < 1; i++){
        ia3 = i*3;
        // for the circular RNA: res[0]-res[n-1] should be excluded from vdW
        if (iCircular){
            k = 1;
        }
        else{
            k = 0;
        }
        // j=i+2 as we exclude two adjacent nucleotides
        // G--A--U--C--'G'--A--C--G--U
        // o--o--o--x--'*'--x--o--o--o
        for (j = i+2; j < natom-k; j++){
            ib3 = j*3;
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > rc2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            r = sqrt(r2);
            invr = 1.0 / r;
            invr2 = invr*invr;
            invr4 = invr2*invr2;
            invr6 = invr4*invr2;
            // accumulate the bond energy and first derivative
            *eng += invr6 * (ljA*invr6 - ljB);
            de = invr4 * invr4 * (fljB - fljA*invr6);
            deriv[ia3] += de*dx;
            deriv[ia3+1] += de*dy;
            deriv[ib3] -= de*dx;
            deriv[ib3+1] -= de*dy;
        }
    }
    */
}

//-----------------------------------------------------------------------------
// eele : computes the van der Waals pairwise interaction using Lennard-Jones
//       potential.
//-----------------------------------------------------------------------------
void eele(prmtop *top, state *sdt, double* E_elec){
    int i, j, k;
    int natom, n, n2;
    int ia, ib, ia3, ib3;
    double dx, dy, de;
    double r2, r, invr2, invr4, invr;
    double *eng = NULL;
    double *deriv = NULL; // first derivative
    double rc2 = top->ecut * top->ecut;
    double charge = top->charge;
    double qA, qB, qq;
    double L_Debye = top->len_Debye; //Debye_length in Angstrom
    //int iCircular = top->iCircular;
    int **verlet_list = sdt->verlet_list;
    int *verlet_natom = sdt->verlet_natom;

    eng = E_elec;
    deriv = sdt->acl;
    natom = sdt->natom;
    qA = charge;
    qB = charge;
    qq = qA*qB;
    for (i = 0; i < natom-2; i++){
        ia3 = i*3;
        ////printf("verlet_natom[%d]=%d\t", i, verlet_natom[i]);
        for (j = 0; j < verlet_natom[i]; j++){
            ////printf("%d  ", verlet_list[i][j]);
            ib3 = verlet_list[i][j] * 3;
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > rc2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            r = sqrt(r2);
            invr = 1.0 / r;
            // accumulate the bond energy and first derivative
            de = qq * invr * exp(-r/L_Debye);
            *eng += de;
            de *= -(1.0/L_Debye + invr) * invr;
            deriv[ia3] += de*dx;
            deriv[ia3+1] += de*dy;
            deriv[ib3] -= de*dx;
            deriv[ib3+1] -= de*dy;
        }
    }
    /*
    // for Excluded list 
    for (i = 0; i < top->nbres-1; i++){
        ia3 = top->nblist[i]*3;
        // j=i+2 as we exclude two adjacent nucleotides
        // G--A--U--C--'G'--A--C--G--U
        // o--o--o--x--'*'--x--o--o--o
        for (j = i+1; j < top->nbres; j++){
            ib3 = top->nblist[j]*3;
            if (ia3-ib3==3 || ib3-ia3==3){
                continue;
            }
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > rc2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            r = sqrt(r2);
            invr = 1.0 / r;
            // accumulate the bond energy and first derivative
            de = qq * invr * exp(-r/L_Debye);
            *eng += de;
            de *= -(1.0/L_Debye + invr) * invr;
            deriv[ia3] += de*dx;
            deriv[ia3+1] += de*dy;
            deriv[ib3] -= de*dx;
            deriv[ib3+1] -= de*dy;
        }
    }
    */
    /*
    // for the circular RNA (not compatible with excluded list yet)
    for (i = 0; i < natom-2; i++){
        ia3 = i*3;
        // for the circular RNA: res[0]-res[n-1] should be excluded from vdW
        if (iCircular){
            k = 1;
        }
        else{
            k = 0;
        }
        // j=i+2 as we exclude two adjacent nucleotides
        // G--A--U--C--'G'--A--C--G--U
        // o--o--o--x--'*'--x--o--o--o
        for (j = i+2; j < natom-k; j++){
            ib3 = j*3;
            // compute displacement
            dx = sdt->crd[ia3] - sdt->crd[ib3];
            dy = sdt->crd[ia3+1] - sdt->crd[ib3+1];
            r2 =  dx*dx + dy*dy;
            if (r2 > rc2){
                continue;   // distance larger than cutoff radius, skip it!
            }
            r = sqrt(r2);
            invr = 1.0 / r;
            // accumulate the bond energy and first derivative
            de = qq * invr * exp(-r/L_Debye);
            *eng += de;
            de *= -(1.0/L_Debye + invr) * invr;
            deriv[ia3] += de*dx;
            deriv[ia3+1] += de*dy;
            deriv[ib3] -= de*dx;
            deriv[ib3+1] -= de*dy;
        }
    }
    */
}

//-----------------------------------------------------------------------------
//
//      5' ---- res1  ----  res2 ---- 3'
//               ||          ||
//      3' ---- res3  ----  res4 ---- 5'        (res: AUGC for 1234)
// 
//  5'-3': AU, UA, GC, CG, GU, UG = 0, 1, 2, 3, 4, 5
//-----------------------------------------------------------------------------
int BpStep(int res1, int res2, int res3, int res4, int iac[]){
    int i5, i3, iBPstep, idx1, idx2, idx3, idx4;
    idx1 = iac[res1];
    //------------------------------
    if (idx1==1){//A
        i5 = 0;//AU
    }
    else if (idx1==4){//C
        i5 = 3;//CG
    }
    else{
        idx3 = iac[res3];
        if (idx1==2){//U
            if (idx3==1){//A
                i5 = 1;//UA
            }
            else{
                i5 = 5;//UG
            }
        }
        else{//G
            if (idx3==2){//U
                i5 = 4;//GU
            }
            else{
                i5 = 2;//GC
            }
        }
    }
    //------------------------------
    idx2 = iac[res2];
    if (idx2==1){//A
        i3 = 0;//AU
    }
    else if (idx2==4){//C
        i3 = 3;//CG
    }
    else{
        idx4 = iac[res4];
        if (idx2==2){//U
            if (idx4==1){//A
                i3 = 1;//UA
            }
            else{
                i3 = 5;//UG
            }
        }
        else{//G
            if (idx4==2){//U
                i3 = 4;//GU
            }
            else{
                i3 = 2;//GC
            }
        }
    }
    //------------------------------
    iBPstep = i5*6+i3;
    return iBPstep;
}
//-----------------------------------------------------------------------------
// Compare: compare the "i,j" with "top, bottom, right, left"
//-----------------------------------------------------------------------------
void Compare(int i, int j, int *topp, int *bottom, int *left, int *right){
    *topp = *topp < i ? *topp : i;
    *bottom = *bottom > i ? *bottom : i;
    *left = *left < j ? *left : j;
    *right = *right > j ? *right : j;
}

//-----------------------------------------------------------------------------
// FindNeighbor: find filled neighboring cmap grids of cmap[i,j]
//-----------------------------------------------------------------------------
void FindNeighbor(int i, int j, int *topp, int *bottom, int *left, int *right,
                int natom, double** cMap, int** NotSearched){
    int max_x, max_y, x, y, n;
    // prevent cross boundary indexing !!
    if (i==0){
        max_x = 2;
    }
    else{
        max_x = 3;
    }
    if (j==natom-1){
        max_y = 2;
    }
    else{
        max_y = 3;
    }
    x = i+1;
    y = j-1;
    for (n = 0; n < max_x; n++){
        if (cMap[x-n][y] && NotSearched[x-n][y]){
            NotSearched[x-n][y] = 0;
            FindNeighbor(x-n, y, topp, bottom, left, right, natom, cMap, NotSearched);
            Compare(x-n, y, topp, bottom, left, right);
        }
    }
    for (n = 1; n < max_y; n++){
        if (cMap[x][y+n] && NotSearched[x][y+n]){
            NotSearched[x][y+n] = 0;
            FindNeighbor(x, y+n, topp, bottom, left, right, natom, cMap, NotSearched);
            Compare(x, y+n, topp, bottom, left, right);
        }
    }
}

//-----------------------------------------------------------------------------
// eBasePair : computes the base-pairing pairwise interaction
//
//            . . . . . . .
//          .               .
//         .                 .
//         .                 .
// (x3,y3) o                 .
//         |                 o (x4,y4)
//         |    Base-Pair    |
// (x2,y2) *--------\        |
//         |         \-------* (x5,y5)
//         |                 |
// (x1,y1) o                 |
//         .                 o (x6,y6)
//         .                 .
//         .                 .
//
//
//-----------------------------------------------------------------------------
void eBasePair(prmtop *top, state *sdt, double* E_basepair, timer *tim, int dpr, FILE *fp, int nstep){
    int i, j, k, l, m, n;
    int natom, n2;
    int ia, ib, ia3, ib3, ntypes=0;
    int idxa, idxb, iparm, iparmm;
    double dx, dy, dz, de;
    double r2, r;
    double *eng = NULL;
    double *deriv = NULL; // first derivative
    double rc2 = top->bcut * top->bcut;
    double eBP, rBP;
    double dx_I, dx_IJ, dx_J;
    double dy_I, dy_IJ, dy_J;
    double rsq_I, rsq_IJ, rsq_J; // rsq means r^2
    int **NotSearched = NULL;
    int iCircular = top->iCircular;
    //
    eng = E_basepair;
    deriv = sdt->acl;
    natom = sdt->natom;
    ntypes = top->pointers[1];
    int *aBPlist5 = NULL, *aBPlist3 = NULL;
    int **BP_partners = NULL, *nBP_partner = NULL;
    int *BPlist5 = NULL, *BPlist3 = NULL, *aStem = NULL;
    int *Stem = NULL, *isTermi = NULL, *BPlistTp = NULL;
    int nBP=0, AnStem=0, nStem=0;
    double **cMap= NULL, **triMap= NULL, **sumMap= NULL;
    int *stlist = NULL, *ndlist = NULL, nBlock=0;
    //
    BP_partners = top->BP_partners;
    nBP_partner = top->nBP_partner;
    //
    aBPlist5 = top->aBPlist5;
    aBPlist3 = top->aBPlist3;
    BPlist5 = top->BPlist5;
    BPlist3 = top->BPlist3;
    BPlistTp = top->BPlistTp;
    //
    isTermi = top->isTermi;
    Stem = top->Stem;
    aStem = top->aStem;
    //
    NotSearched = sdt->NotSearched;
    cMap = sdt->cMap;
    triMap = sdt->triMap;
    sumMap = sdt->sumMap;
    stlist = sdt->stlist;
    ndlist = sdt->ndlist;
    
    int minLoop = top->min_loop;
    int minbpFLAG = top->minbpFLAG;
    int minStembp = top->minStembp;
    double bulge_p = top->bulge_p;

    //
    for (i = 0; i < natom; i++){
        isTermi[i] = 0;
        aStem[i] = 1;
        Stem[i] = 1;
        stlist[i] = 0;
        ndlist[i] = 0;
        ///
        BP_partners[0][i] = 0;
        BP_partners[1][i] = 0;
        nBP_partner[i] = 0;
        /*
        for (j = i; j < natom; j++){
            NotSearched[i][j] = 1;
            NotSearched[j][i] = 1;
            cMap[i][j] = 0.0;       // necessary!
            cMap[j][i] = 0.0;       // necessary!
            //triMap[i][j] = 0.0;     // maybe not necessary
            //sumMap[i][j] = 0.0;     // maybe not necessary
        }
        */
    }
tim->t_BP_others += punchClock(tim->tt);

    //--------------------------------------------------------------------------------------
    //------ Step 1. Calculate the contact map ---------------------------------------------
    //--------------------------------------------------------------------------------------
    //// if dpr == 2, it's constrained-MD where only near-NS cMap values require calculation
    //printf("dpr=%d\n", dpr);
    if (dpr == 2){
        int *dprBplist=NULL, dpr_nBP;
        dprBplist = top->dprBplist;
        dpr_nBP = top->dpr_nBP;
        //printf("dpr_nBP=%d\n", dpr_nBP);
        for (n = 0; n < dpr_nBP; n++){
            // get atom index and force field parm index
            // i-j is res5-res3 of a basepair
            i = dprBplist[n*2];
            // ib could be [j-2],[j-1], [j], [j+1], and [j+2]
            for (m = -2; m < 3; m++){
                j = dprBplist[n*2+1] + m;
                if (j<0 || j>natom-1){
                    continue;
                }
                ia3 = i*3;
                ib3 = j*3;
                //
                idxa = top->iac[i]-1;
                idxb = top->iac[j]-1;
                iparm = top->ico[ntypes*idxa + idxb] -1;
                //
                dx_IJ = sdt->crd[ib3] - sdt->crd[ia3];
                dy_IJ = sdt->crd[ib3+1] - sdt->crd[ia3+1];
                rsq_IJ =  dx_IJ*dx_IJ + dy_IJ*dy_IJ;
                if (rsq_IJ > rc2){
                    continue;   // distance > cutoff, skip it!
                }
                //  not Sqrt() yet! For ease of calculation
                cMap[i][j] = top->bpFLAG[iparm] / rsq_IJ;   // considering pairing type yet!
                cMap[j][i] = cMap[i][j];
                //printf("ia=%d, ib=%d, rsq_IJ=%lf, cMap[%d][%d]=%lf, bpFLAG=%d\n", i, j, rsq_IJ, i,j,top->bpFLAG[iparm] / rsq_IJ, top->bpFLAG[iparm]);
            }
        }
    }
    /*
    else{// normal MD
        for (i = 0; i < natom-minLoop-1; i++){  // tri-loop
            ia3 = i*3;
            idxa = top->iac[i]-1;
            // ------------- 0 --- 1 --- 2 --- 3 --- 4 --- 5 --- 6 -----
            // vdw:          i     -    [o     o     o     o     o
            // BP(min=3):    i     -     x     x    [o     o     o
            // BP(min=4):    i     -     x     x     x    [o     o
            // j = minLoop-1: because the 'x' above should be eliminated from verlet list
            for (j = minLoop-1; j < verlet_natom[i]; j++){
                ib = verlet_list[i][j];
                ib3 = ib*3;
                idxb = top->iac[ib]-1;
                iparm = top->ico[ntypes*idxa + idxb] -1;
                // check whether this pair is BP-forming 'AU'/'GC'/'GU'
                ////if (top->bpFLAG[iparm]==0){
                ////    continue;
                ////}
                // compute displacement between two BP-forming beads
                dx_IJ = sdt->crd[ib3] - sdt->crd[ia3];
                dy_IJ = sdt->crd[ib3+1] - sdt->crd[ia3+1];
                rsq_IJ =  dx_IJ*dx_IJ + dy_IJ*dy_IJ;
                if (rsq_IJ > rc2){
                    continue;   // distance > cutoff, skip it!
                }
                //  not Sqrt() yet! For ease of calculation
                cMap[i][ib] = top->bpFLAG[iparm] / rsq_IJ;   // considering pairing type yet!
                cMap[ib][i] = cMap[i][ib];
                //printf("cMap[%d][%d]=%lf\t", i, j, cMap[i][j]);
            }
        }
    }
    */

    // For the circular RNA: zero some cmap values related to termini
    if (iCircular){
        for (i = 0; i < minLoop; i++){
            for (j = 0; j < minLoop-i; j++){
                cMap[i][natom-1-j] = 0.0;
            }
        }
    }
//tim->t_BP_DP += punchClock(tim->tt);

    /*
    //----------------------------------------------------------------------------------------------------
    // ------ Step 2. Split into the blocks and find out the bolck harboring the maximum cmap ------------
    //----------------------------------------------------------------------------------------------------
    int *max_iblock = sdt->max_iblock;
    int **col_iblock = sdt->col_iblock;
    int **row_iblock = sdt->row_iblock;
    double max_row;
    int iblock, isBreak;
    for (n = 0; n < 2; n++){
        for (i = 0; i < natom; i++){
            max_row = 0.0;
            iblock = 1;
            isBreak = 0;
            //-------- vetical: mirrored cmap ------------
            for (j = 0; j < i-1; j++){
                if (cMap[j][i]){
                    if (cMap[j+1][i] < 0.00001){
                        isBreak = 1;
                    }
                    else{
                        isBreak = 0;
                    }
                    //
                    row_iblock[i][j] = iblock;
                    col_iblock[j][i] = iblock;
                    //
                    if (max_row < cMap[j][i]){
                        max_row = cMap[j][i];
                        max_iblock[i] = row_iblock[i][j];
                    }
                    if (isBreak){
                        iblock += 1;
                    }
                }
            }
            //-------- parallel: valued cmap ------------
            for (j = i; j < natom-1; j++){
                if (cMap[i][j]){
                    if (cMap[i][j+1] < 0.00001){
                        isBreak = 1;
                    }
                    else{
                        isBreak = 0;
                    }
                    //
                    row_iblock[i][j] = iblock;
                    col_iblock[j][i] = iblock;
                    //
                    if (max_row < cMap[i][j]){
                        max_row = cMap[i][j];
                        max_iblock[i] = row_iblock[i][j];
                    }
                    if (isBreak){
                        iblock += 1;
                    }
                }
            }
            // the last column
            if (cMap[i][natom-1]){
                row_iblock[i][natom-1] = iblock;
                col_iblock[natom-1][i] = iblock;
                if (max_row < cMap[i][natom-1]){
                    max_row = cMap[i][natom-1];
                    max_iblock[i] = row_iblock[i][natom-1];
                }
            }
            //printf("max_row[%3d]=%6.4lf, max_iblock[%3d]=%2d\n", i, max_row, i, max_iblock[i]);
        }
        //test
        //---------------------------------------------------------------
        for (i = 0; i < natom; i++){
            for (j = i; j < natom; j++){
                if (cMap[i][j]){
                    printf("row_iblock[%d][%d]=%d  ", i, j, row_iblock[i][j]);
                    printf("col_iblock[%d][%d]=%d  ", i, j, col_iblock[i][j]);
                    printf("max_iblock[%d]=%d, max_iblock[%d]=%d\n", i, max_iblock[i], j, max_iblock[j]);
                }
            }
            printf("\n");
        }
        
        printf("Before triple-bp exclusion\n");
        for (i = 0; i < natom; i++){
            for (j = i; j < natom; j++){
                if (cMap[i][j]){
                    printf("cMap[%d][%d]=%.4lf  ", i, j, cMap[i][j]);
                }
            }
            printf("\n");
        }
        //---------------------------------------------------------------
    
    
        //test
        
        //----------------------------------------------------------------------------------------------------
        //------ Step 3. Prevention of Triple-BP (using only half of the cMap) --------------------------
        // This step also eliminates the PK-forming blocks, by zeroing non-maximum cmap values
        //----------------------------------------------------------------------------------------------------
        for (i = 0; i < natom; i++){
            for (j = i; j < natom; j++){
                if (cMap[i][j]){
                    // check whether the belonging block of this cmap is the maximum of row AND column !
                    if (row_iblock[i][j]!=max_iblock[i] || col_iblock[i][j]!=max_iblock[j]){
                        cMap[i][j] = 0.0;
                    }
                }
            }
        }
        //test
        //---------------------------------------------------------------
        printf("After triple-bp exclusion\n");
        for (i = 0; i < natom; i++){
            for (j = i; j < natom; j++){
                if (cMap[i][j]){
                    printf("cMap[%d][%d]=%.4lf  ", i, j, cMap[i][j]);
                }
            }
            printf("\n");
        }
        //---------------------------------------------------------------
    }
tim->t_BP_exTri += punchClock(tim->tt);
    */
    
    
    //test

    //----------------------------------------------------------------------------------------------------
    //------ Step 2. Select the stem blocks according to contact map (recursively) -------------------------
    //----------------------------------------------------------------------------------------------------
    int topp, bottom, right, left;  // 'topp' not to collide with '*top'
    nBlock = 0;
    for (i = 0; i < natom-3; i++){
        for (j = natom-1; j > i+2; j--){
            //printf("cMap[%d][%d]=%lf\t", i, j, cMap[i][j]);
            if (cMap[i][j] && NotSearched[i][j]){
                topp = i;
                bottom = i;
                right = j;
                left = j;
                NotSearched[i][j] = 0;
                FindNeighbor(i, j, &topp, &bottom, &left, &right, natom, cMap, NotSearched);
                ///printf("i=%d, j=%d, m=%d, n=%d\n", i, j, m, n);
                // cull out stem block smaller than 2x2, nammely 1x1, 1xN, Nx1
                if ((bottom-topp) && (right-left)){
                    ndlist[2*nBlock] = topp;
                    ndlist[2*nBlock+1] = right;
                    stlist[2*nBlock] = bottom;
                    stlist[2*nBlock+1] = left;
                    nBlock += 1; // the number of stem-blocks
                }
            }
        }
    }
    //printf("nBlock=%d\n", nBlock);

    //----------------------------------------------------------------------------------------------------
    //------ Step 3. Dynamic Progromming to find best-Stem of each stem-blosk ------------------------
    //----------------------------------------------------------------------------------------------------
    int i0, j0, it, jt, LEN, iA, jA, ci, cj, max_i, max_j;
    double csum, max_BP;
    nBP = 0;
    for (m = 0; m < nBlock; m++){
        i0 = ndlist[2*m];
        j0 = ndlist[2*m+1];
        it = stlist[2*m];
        jt = stlist[2*m+1];
        LEN = it - i0 + j0 - jt + 1;  //'+1' is important!!!
        // In case of continuous hairpin! that a helix traverse the whole cmap
	    // and the Trimap/Summap will exceed the maximum size
	    LEN = LEN < natom-2 ? LEN : natom-2;
	    iA = i0;
        jA = j0 - LEN +1;   //'+1' is important!!!
        ////printf("\n\nLEN=%d\n\n", LEN);
        ////printf("i0=%d, j0=%d, it=%d, jt=%d, iA=%d, jA=%d\n\n", i0, j0, it, jt, iA, jA);
        // we fill Map[LEN][LEN] to the Top-Left triangle of 'triMap'
        for (i = 0; i < LEN+2; i++){
            for (j = 0; j < i+1; j++){
                triMap[i][j] = 0.0;
                sumMap[i][j] = 0.0;
            }
        }
        for (i = 0; i < LEN; i++){
            n = 0;
            for (j = LEN-i-1; j < LEN; j++){
                triMap[i][n] = cMap[n+iA][j+jA];
                n += 1;
            }
        }
        //test
	/*
        printf("\n\n Trimap \n\n");
        for (i = 0; i < LEN; i++){
            for (j = 0; j < LEN; j++){
                printf("triMap[%d][%d]=%.4lf  ", i, j, triMap[i][j]);
            }
            printf("\n");
        }
        */
        //test
        //--------------- **** Summation  **** -------------------------------------
        n = LEN+2;
        for (i = n-3; i > -1; i--){
            for (j = 0; j < i+1; j++){
                sumMap[i][j] = sumMap[i+1][j] > sumMap[i+1][j+1] ?
                                (sumMap[i+1][j]- bulge_p) : (sumMap[i+1][j+1]- bulge_p);
                // bulge penalty !!
                //sumMap[i][j] -= 0.005;
                //
                sumMap[i][j] = sumMap[i+2][j+1]+triMap[i][j] > sumMap[i][j] ? 
                                sumMap[i+2][j+1]+triMap[i][j] : sumMap[i][j];
            }
            ////printf("%d, %.4lf\n", i, sumMap[i][j-1]);
        } 
        ////printf("sumMap[0][0]=%lf\n", sumMap[0][0]);
        //--------------- **** Where-to-start ? **** -------------------------------------
        max_i = 0;
        max_j = 0;
        max_BP = 0.0;
        for (i = 0; i < LEN; i++){
            if (sumMap[i][0]>max_BP){
                max_BP = sumMap[i][0];
                max_i = i;
            }
            if (sumMap[i][i]>max_BP){
                max_BP = sumMap[i][i];
                max_i = i;
                max_j = i;
            }
        }
        ////printf("sumMap[-0][-0]=%lf, max_BP=%lf\n", sumMap[0][0], max_BP);
        ////printf("max_i=%d, max_j=%d\n", max_i, max_j);
        //------------------ Backtracking ----------------------------------------------
        ci = max_i;
        cj = max_j;
        while (ci < LEN){
            csum = sumMap[ci][cj];
            if (csum == (sumMap[ci+1][cj] - bulge_p)){
                ci += 1;
            }
            else if (csum == (sumMap[ci+1][cj+1] - bulge_p)){
                ci += 1;
                cj += 1;
            }
            else{// csum == sumMap[ci+2][cj+1]+triMap[ci][cj]
                // The stem triangle block starts at A(xA,yA)
                // '-2' for 'j' is important, because we add two
                // additional row not belonging to original ctMap
                // but it cannot be added for 'i'!!!
                ia = cj + iA;
                ib = n-1 - ci + cj + jA - 2;
                if (abs(ia-ib) > minLoop){
                    idxa = top->iac[ia]-1;
                    idxb = top->iac[ib]-1;
                    iparm = top->ico[ntypes*idxa + idxb] -1;
                    ////printf("%d, %d, %d\n", ia, ib, iparm);
                    if (top->bpFLAG[iparm]>minbpFLAG && cMap[ia][ib]>0.0){  // AU, GU, GC, or GA
		    // 'cMap[ia][ib]>0.0' is for the cmap zeroed during triple-Bp exclusion,
		    // which should no longer be considerer in this stem block
                        aBPlist5[nBP] = ia;
                        aBPlist3[nBP] = ib;
                        nBP += 1;
			////printf("%d, %d, %d\n", ia, ib, iparm);			
                    }
                }
                ci += 2;
                cj += 1;
            }
        }
    }
    /*
    for (i = 0; i < nBP; i++){
        printf("%d-%d\t", aBPlist5[i]+1, aBPlist3[i]+1);
    }
    */
    //printf("nBP=%d\n\n", nBP);
    

    //----------------------------------------------------------------------------------------------------
    //------------- Step 3.2. trim the BP list into stem segments --------------------------
    //----------------------------------------------------------------------------------------------------
    int next;
    AnStem = 0;
    next = 1;
    for (i = 0; i < nBP-1; i++){
        if (aBPlist5[i+1]-aBPlist5[i] == 1 &&
            aBPlist3[i]-aBPlist3[i+1] == 1){
            aStem[AnStem] += 1;
            next = 0;
        }
        else{
            next = 1;
        }
        if (next){
            AnStem += 1;
        }
    }
    if (nBP){
        AnStem += 1;
    }
tim->t_BP_DP += punchClock(tim->tt);


    //----------------------------------------------------------------------------------------------------
    //------------- Step 4. Exclude the parallel stem, trans-isomer, and Pseudo-Knot ---------------------
    //----------------------------- Stem Exclusion: Vector Method --------------------------------------
    //----------------------------------------------------------------------------------------------------
    int *orient=NULL, Accept;
    double vecAB[2], vecCD[2], dot_ABCD;    // for parallel-stem exclusion use
    double moveX, moveY, crdTx, crdTy, crdMx, crdMy, theta, crdMrotY;
    int upward, downward;
    upward = 0;
    downward = 0;
    orient = top->orient;
    n = 0;  // index for aBPist
    //-------------------------------------------------------------------------------------
    //------------- Step 4.1. Check parallel-stem  ----------------------------------------
    //-------------------------------------------------------------------------------------
    /// whether there is any PK?
    for (i = 0; i < AnStem; i++){
        // -------  consider all the basepairs including lonepairs --------
        // check the unwanted "parallel" 2bp-stem  or lone-pair
        if (aStem[i]<3){    // 2-bp stem or lone-pair
            //-----------------------------------------------------------------------------
            // parallel stem can only occur in 2-bp stem as the restriction of bp-cutoff
            // only need the first basepair, as lonepair has no second basepair
            //
            //           _---_
            //          |     |
            //          *     *
            //          |     |
            //          B === C  <--- (second basepair, optional)
            //          |     |
            //          A === D  <-- (first basepair)
            //          |     |
            //      5'--*     *--3'
            //
            //-----------------------------------------------------------------------------
            for (j = 0; j < 2; j++){    // x and y
                vecAB[j] = sdt->crd[3* (aBPlist5[n]+1) +j] - sdt->crd[3 * aBPlist5[n] +j];
                vecCD[j] = sdt->crd[3* aBPlist3[n] +j] - sdt->crd[3 * (aBPlist3[n]-1) +j];
            }
            // the dot of vector_AB and vector_CD
            dot_ABCD = vecAB[0]*vecCD[0] + vecAB[1]*vecCD[1];
            if (dot_ABCD > 0.0){    // angle(AB,CD) < 90 degreed, parallel
                orient[i] = 0;
                n += aStem[i];
                // No need of further judging it's trans or PK, just exclude it right now
                continue;
            }
        }
        //---------------------------------------------------------------------
        //--------- Step 4.2. Check the trans-isomer and Pseudoknots -----------
        //---------------------------------------------------------------------
        // 1. translation to superimpose 'bp5-5' onto the origin
        moveX = sdt->crd[3*aBPlist5[n]];
        moveY = sdt->crd[3*aBPlist5[n]+1];
        // translate the crd of "bp5-3" 
        crdTx = sdt->crd[3*aBPlist3[n]] - moveX;
        crdTy = sdt->crd[3*aBPlist3[n]+1] - moveY;
        // get and translate the middle of "second basepair"
        // NOTICE: 'second-bp' is replace by the adjacent [i+1] and [j-i]
        //          for the harmony of stems and lonepairs (without 2nd bp)
        crdMx = (sdt->crd[3* (aBPlist5[n]+1)   ] + sdt->crd[3* (aBPlist3[n]-1)   ])/2 - moveX;
        crdMy = (sdt->crd[3* (aBPlist5[n]+1) +1] + sdt->crd[3* (aBPlist3[n]-1) +1])/2 - moveY;
        // 2. rotation 'bp5' to align with "x+" axis
        theta = acos(crdTx / sqrt(crdTx*crdTx + crdTy*crdTy));  //rotation angle
        if (crdTy>0){   //change to clockwise rotation!
            theta = -theta;
        }
        // Rotation Matrix = [[cos(theta), -sin(theta)],[sin(theta), cos(theta)]]
        crdMrotY = crdMx*sin(theta) + crdMy*cos(theta);
        crdMrotY = (crdMrotY<-1.0) ? -1.0 : (crdMrotY>1.0) ? 1.0 : crdMrotY;
        // 3. determine the orientation of stem
        if (crdMrotY>0.0){ // upwards
            orient[i] = 1;
            upward += aStem[i];
        }
        else{       // downwards
            orient[i] = -1;
            downward += aStem[i];
        }
        n += aStem[i];
    }
    //----------------------------------------------------------------------
    //--------- Step 4.3. Compare and determine the more stable orientation
    //---------            (the base pairs are not adopted yet at this step)
    //----------------------------------------------------------------------
    if (upward>downward){
        Accept = 1;
    }
    else{
        Accept = -1;
    }
    // Do the acception:
    n = 0;
    for (i = 0; i < AnStem; i++){
        //printf("aStem[%d]=%d, orient[%d]=%d\n", i, aStem[i], i, orient[i]);
        if (orient[i] == Accept){
            // abandon the lone-pairs
            if (aStem[i] < minStembp){
                n += 1;
                continue;
            }
            // finalize the basepairs
            for (j = 0; j < aStem[i]; j++){
                //printf("BP-%d\n", n);
                //BPlist5[nBP] = aBPlist5[n];
                //BPlist3[nBP] = aBPlist3[n];
                //-------the ct list (for triplet exclusion)--------------
                if (nBP_partner[aBPlist5[n]] > 1){// quadruplet
                    fprintf(fp, "Warning! quadruplet at frame-%d, res-%d\t", nstep, aBPlist5[n]);
                    fprintf(fp, "(%d-%d)\n", aBPlist5[n]+1, aBPlist3[n]+1);
                    n += 1;
                    continue;
                }
                if (nBP_partner[aBPlist3[n]] > 1){// quadruplet
                    fprintf(fp, "Warning! quadruplet at frame-%d, res-%d\t", nstep, aBPlist3[n]);
                    fprintf(fp, "(%d-%d)\n", aBPlist5[n]+1, aBPlist3[n]+1);
                    n += 1;
                    continue;
                }
                //-----------------------------------------------------------------
                k = nBP_partner[aBPlist5[n]];
                ///printf("k=%d, aBPlist5[%d]=%d\n", k, n, aBPlist5[n]);
                BP_partners[k][aBPlist5[n]] = aBPlist3[n] + 1;
                nBP_partner[aBPlist5[n]] += 1;
                ///printf("BP_partners[%d][%d]=%d\n", k, aBPlist5[n], aBPlist3[n] + 1);
                ///printf("nBP_partner[%d]=%d\n", aBPlist5[n], nBP_partner[aBPlist5[n]]);
                //
                k = nBP_partner[aBPlist3[n]];
                ///printf("k=%d, aBPlist3[%d]=%d\n", k, n, aBPlist3[n]);
                BP_partners[k][aBPlist3[n]] = aBPlist5[n] + 1;
                ///printf("BP_partners[%d][%d]=%d\n", k, aBPlist3[n], aBPlist5[n] + 1);
                nBP_partner[aBPlist3[n]] += 1;
                ///printf("nBP_partner[%d]=%d\n", aBPlist3[n], nBP_partner[aBPlist3[n]]);
                //--------------------------------------------------------
                n += 1;
            }
        }
        // abandon lonepair and stems from opposite orientation
        else{
            n += aStem[i];
        }
    }
    //top->nStem = nStem;
    //top->nBP = nBP;
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
tim->t_BP_exPK += punchClock(tim->tt);


    //test
    /*
    printf("\nPseudoknot trimming!\n");
    printf("nBP=%d\n", nBP);
    for (i = 0; i < nBP; i++){
        printf("%d-%d\t", aBPlist5[i]+1, aBPlist3[i]+1);
    }
    printf("AnStem=%d\n", AnStem);
    printf("\n upward=%lf, downward=%lf\n", upward, downward);
    */
    //-----------------------------------------------------------------------------
    //------------- Step 5.1. Exclude the triplets ----------------------------------
    //-----------------------------------------------------------------------------
    /// check
    double dist_1, dist_2;
    /*
    for (i = 0; i < natom; i++){
        printf("nBP_partner[%d]=%d, BP_partners[0][%d]=%d, BP_partners[1][%d]=%d\n",
            i, nBP_partner[i], i, BP_partners[0][i], i, BP_partners[1][i]);
    }
    */
    for (i = 0; i < natom; i++){
        if (nBP_partner[i] > 1){
            ia = BP_partners[0][i]-1;
            ib = BP_partners[1][i]-1;
            dist_1 = cMap[i][ia];
            dist_2 = cMap[i][ib];
            //
            nBP_partner[i] -= 1;
            if (dist_1 > dist_2){
                // keep the fist bp, and vacate the second bp
                BP_partners[1][i] = 0;
                // rearange the ib
                BP_partners[0][ib] = BP_partners[1][ib];
                nBP_partner[ib] -= 1;
                // (if ib has only one bp, BP_partners[1][ib]=0,
                // if ib has two bps, BP_partners[0][ib] must be a smaller index than ib
                // and the we should replace BP_partners[0][ib] with BP_partners[1][ib])
            }
            else{
                // replace the first with the second, then vacate the second bp
                BP_partners[0][i] = BP_partners[1][i];
                BP_partners[1][i] = 0;
                // rearange the ia,
                BP_partners[0][ia] = BP_partners[1][ia];
                nBP_partner[ia] -= 1;
                // (if ia has only one bp, BP_partners[1][ia]=0,
                // if ia has two bps, BP_partners[0][ia] must be a smaller index than ia
                // and the we should replace BP_partners[0][ia] with BP_partners[1][ia])
            }
            
        }
    }

    //----------------------------------------------------------------------------------------------------
    //------------- Step 5.2. trim the BP list into stem segments ----------------------------------------
    //----------------------------------------------------------------------------------------------------
    // ---- Re-select the base-pairs ------------------------
    //int res55, res53;
    //int res35, res33;
    nBP = 0;
    for (i = 0; i < natom; i++){
        if (nBP_partner[i]){
            BPlist5[nBP] = i;
            BPlist3[nBP] = BP_partners[0][i] - 1;
            // prevent double selection
            nBP_partner[BP_partners[0][i]-1] -= 1;
            nBP += 1;
        }
    }
    top->nBP = nBP;
    /*
    printf("There are %d BasePairs:   \n", nBP);
    for (i = 0; i < nBP; i++){
        printf("%d: %d-%d (%d)\n", i, BPlist5[i], BPlist3[i], isTermi[i]);
    }
    */
    

    //------ Re-trim the stems ------------------------------
    next = 1;
    nStem = 0;
    for (i = 0; i < nBP-1; i++){
        if (BPlist5[i+1]-BPlist5[i] == 1 &&
            BPlist3[i]-BPlist3[i+1] == 1){
            Stem[nStem] += 1;
            next = 0;
        }
        else{
            next = 1;
        }
        if (next){
            nStem += 1;
        }
    }
    if (nBP){
        nStem += 1;
    }
    top->nStem = nStem;

    /*
    printf("There are %d Stems:\n", nStem);
    for (i = 0; i < nStem; i++){
        printf("  length of Stem %d is %d\n", i, Stem[i]);
    }
    */
tim->t_BP_exTri += punchClock(tim->tt);

    //-----------------------------------------------------------------------------
    //------------- Step 5.3. Mark the Termini of each stem -----------------------
    //-----------------------------------------------------------------------------
    n = 0;
    k = 0;
    for (i = 0; i < nStem; i++){
        // re-check of lone-pairs (results from some 2bp stems before triplet-exclusion)
        m = Stem[i];
        if (m+1 > minStembp){
            for (j = 0; j < m; j++){
                aBPlist5[k] = BPlist5[n+j];
                aBPlist3[k] = BPlist3[n+j];
                k += 1;
            }
        }
        isTermi[n] = 5;
        isTermi[n+m-1] = 3;
        n += m;
    }
    top->nBP = k;
    // Note!!! top->nBP isn't nBP, the former excluded the lonepair(if minStembp=2)


    //-----------------------------------------------------------------------------
    //------------- Step 9. Determine the stacking and angle list -----------------
    //-----------------------------------------------------------------------------
    int a, b, c, d;
    int idx, *anglist = NULL, *stklist = NULL;
    int isStk, nangle, maxNstk, nSTK, iBPstep;
    maxNstk = 4*4*natom;
    anglist = top->anglist;
    stklist = top->stklist;
    for (i = 0; i < maxNstk; i++){
        stklist[i] = 0;
    }
    for (i = 0; i < 2*natom; i++){
        anglist[i] = 0;
    }
    top->nangle = 0;
    nangle = 0;
    nSTK = 0;
    isStk = 0;
    top->nSTk = 0;
    for (i = 0; i < nBP; i++){
        if (isTermi[i]==3){
            // this code deals perfectly with lonepairs, which are always marked as 3'-stem
            isStk = 0;
        }
        if (isStk==1){
            anglist[nangle] = BPlist5[i] * 3;
            anglist[nangle+1] = BPlist3[i] * 3;
            nangle += 2;
            //
            iBPstep = BpStep(BPlist5[i], BPlist5[i+1], BPlist3[i], BPlist3[i+1], top->iac);
            idx = 4*nSTK;
            a = BPlist3[i] * 3;
            b = BPlist5[i] * 3;
            c = BPlist3[i+1] * 3;
            d = BPlist5[i+1] * 3;
            //printf("a=%d, b=%d, c=%d, d=%d\n", a, b, c, d);
            stklist[idx] = a;     // atom-i
            stklist[idx+1] = b;     // atom-j (center)
            stklist[idx+2] = d;     // stom-k
            stklist[idx+3] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+4] = b;     // atom-i
            stklist[idx+5] = a;     // atom-j (center)
            stklist[idx+6] = c;     // stom-k
            stklist[idx+7] = iBPstep;
            stklist[idx+8] = b;     // atom-i
            stklist[idx+9] = d;     // atom-j (center)
            stklist[idx+10] = c;     // stom-k
            stklist[idx+11] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+12] = a;     // atom-i
            stklist[idx+13] = c;     // atom-j (center)
            stklist[idx+14] = d;     // stom-k
            stklist[idx+15] = iBPstep;
            //
            nSTK += 4;
        }
        // this code deals perfectly with lonepairs, which are always marked as isTermi[n]=3
        if (isTermi[i]==5){
            isStk = 1;
            iBPstep = BpStep(BPlist5[i], BPlist5[i+1], BPlist3[i], BPlist3[i+1], top->iac);
            idx = 4*nSTK;
            a = BPlist3[i] * 3;
            b = BPlist5[i] * 3;
            c = BPlist3[i+1] * 3;
            d = BPlist5[i+1] * 3;
            //printf("a=%d, b=%d, c=%d, d=%d\n", a, b, c, d);
            stklist[idx] = a;     // atom-i
            stklist[idx+1] = b;     // atom-j (center)
            stklist[idx+2] = d;     // stom-k
            stklist[idx+3] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+4] = b;     // atom-i
            stklist[idx+5] = a;     // atom-j (center)
            stklist[idx+6] = c;     // stom-k
            stklist[idx+7] = iBPstep;
            stklist[idx+8] = b;     // atom-i
            stklist[idx+9] = d;     // atom-j (center)
            stklist[idx+10] = c;     // stom-k
            stklist[idx+11] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+12] = a;     // atom-i
            stklist[idx+13] = c;     // atom-j (center)
            stklist[idx+14] = d;     // stom-k
            stklist[idx+15] = iBPstep;
            //
            nSTK += 4;
        }
    }
    top->nangle = nangle;
    top->nSTk = nSTK;

    // test
    /*
    printf("There are %d Stems:\n", nStem);
    for (i = 0; i < nStem; i++){
        printf("  length of Stem %d is %d\n", i, Stem[i]);
    }
    printf("There are %d BasePairs:   \n", nBP);
    for (i = 0; i < nBP; i++){
        printf("%d: %d-%d (%d)\n", i, BPlist5[i], BPlist3[i], isTermi[i]);
    }
    printf("\n");
    printf("There are %d Angles:\n", top->nangle);
    for (i = 0; i < top->nangle; i++){
        printf("%d\n", anglist[i]);
    }
    printf("There are %d Stacks:\n", top->nSTk);
    for (i = 0; i < top->nSTk; i++){
        printf("%d: %d-%d-%d\n", i, stklist[i*4], stklist[i*4+1], stklist[i*4+2]);
    }
    */
    
tim->t_BP_Angs += punchClock(tim->tt);


    //------------- mark 5' and 3' termini of the RNA strand ---------------------
    // Notice! can't mark this before stacking-determination, because 5' and 3' BP can be lone-pair
    // this is only for "index-out-of-range" consideration when calculating BP-theta angle
    /*
    for (i = 0; i < nBP; i++){
        if (BPlist5[i] == 0 || BPlist3[i] == natom-1){
            isTermi[i] = 5;
        }
    }
    */

    //------------- Step 10. Apply harmonic restraint to basepairs ------------------
    double bphm_eq, bphm_k, bphm_p, dl;
    int nHB;
    bphm_eq = top->bphm_eq;
    bphm_k = top->bphm_k;
    bphm_p = top->bphm_p;
    n = 0;
    //printf("Restraint on: ");
    for (i = 0; i < top->nBP; i++){
        // get atom index and force field parm index
        //printf("%d-%d ", BPlist5[i], BPlist3[i]);
        ia = aBPlist5[i]*3;
        ib = aBPlist3[i]*3;
        // get the number of HB for different BP
        idxa = top->iac[aBPlist5[i]]-1;
        idxb = top->iac[aBPlist3[i]]-1;
        iparm = top->ico[ntypes*idxa + idxb] -1;
        nHB = top->nHB[iparm];
        // Compute displacement
        dx = sdt->crd[ia] - sdt->crd[ib];
        dy = sdt->crd[ia+1] - sdt->crd[ib+1];
        //
        r2 =  dx*dx + dy*dy;
        r = sqrt(r2);
        dl = r - bphm_eq;
        //printf("i=%d, ia=%d, ib=%d, dl=%lf\n", i, ia, ib, dl);

        // Accumulate the bond energy and first derivative                       // Notice: it's gradient, not force yet
        *eng += nHB * (bphm_k * dl * dl + bphm_p);
        de = 2.0 * nHB * bphm_k * dl / r;
        deriv[ia] += de*dx;
        deriv[ia+1] += de*dy;
        deriv[ib] -= de*dx;
        deriv[ib+1] -= de*dy;
    }

    //printf("\n");

tim->t_BP_others += punchClock(tim->tt);
    // ------------------------------------------------------------
}

//-----------------------------------------------------------------------------
// eBasePair_noDP : computes the base-pairing pairwise interaction without
//                  using Dynamic Programming.
//                  Note: Requiring 2D structure prescribed.
//                        Only exact basepairs are checked (within the bcut)
//-----------------------------------------------------------------------------
void eBasePair_noDP(prmtop *top, state *sdt, double* E_basepair, timer *tim){
    int i, j, k, m, n, natom;
    int ia, ib, iparm;
    double rc2 = top->bcut * top->bcut;
    double dx, dy, de;
    double r2, r, dl;
    double *eng = NULL;
    double *deriv = NULL;
    //
    eng = E_basepair;
    deriv = sdt->acl;
    natom = sdt->natom;

    

    //----------------------------------------------------------------------------------------------------
    //------------- Step 1. trim the BP list into stem segments ----------------------------------------
    //----------------------------------------------------------------------------------------------------
    // ---- Select the base-pairs from prescribed 2D structure ------------------------
    int *dprBplist=NULL, nBP = 0;
    int *BPlist5 = NULL, *BPlist3 = NULL;
    int *aBPlist5 = NULL, *aBPlist3 = NULL;
    BPlist5 = top->BPlist5;
    BPlist3 = top->BPlist3;
    aBPlist5 = top->aBPlist5;
    aBPlist3 = top->aBPlist3;
    dprBplist = top->dprBplist;
    /*
    printf("There are confined %d BasePairs:   \n", top->dpr_nBP);
    for (i = 0; i < top->dpr_nBP; i++){
        printf("%d: %d-%d\t", i, rBplist[i*2], rBplist[i*2+1]);
    }
    printf("\n");
    */
    nBP = 0;
    for (i = 0; i < top->dpr_nBP; i++){
        ia = dprBplist[i*2];
        ib = dprBplist[i*2+1];
        // Compute displacement
        dx = sdt->crd[ia*3] - sdt->crd[ib*3];
        dy = sdt->crd[ia*3+1] - sdt->crd[ib*3+1];
        r2 =  dx*dx + dy*dy;
        // check whether the basepairs are within the cutoff
        if (r2 > rc2){
            continue;
        }
        BPlist5[nBP] = ia;
        BPlist3[nBP] = ib;
        nBP += 1;
    }
    top->nBP = nBP;
    /*
    printf("There are %d BasePairs:   \n", nBP);
    for (i = 0; i < nBP; i++){
        printf("%d: %d-%d\t", i, BPlist5[i], BPlist3[i]);
    }
    printf("\n");
    */

    //------ Re-trim the stems ------------------------------
    int next, nStem = 0;
    int *Stem = NULL;
    int *isTermi = NULL;
    isTermi = top->isTermi;
    Stem = top->Stem;
    for (i = 0; i < natom; i++){
        isTermi[i] = 0;
        Stem[i] = 1;
    }
    for (i = 0; i < nBP-1; i++){
        if (BPlist5[i+1]-BPlist5[i] == 1 &&
            BPlist3[i]-BPlist3[i+1] == 1){
            Stem[nStem] += 1;
            next = 0;
        }
        else{
            next = 1;
        }
        if (next){
            nStem += 1;
        }
    }
    if (nBP){
        nStem += 1;
    }
    top->nStem = nStem;

    //-----------------------------------------------------------------------------
    //------------- Step 2. Mark the Termini of each stem -----------------------
    //-----------------------------------------------------------------------------
    int minStembp = top->minStembp;
    n = 0;
    k = 0;
    for (i = 0; i < nStem; i++){
        // re-check of lone-pairs (results from some 2bp stems before triplet-exclusion)
        m = Stem[i];
        if (m+1 > minStembp){
            for (j = 0; j < m; j++){
                aBPlist5[k] = BPlist5[n+j];
                aBPlist3[k] = BPlist3[n+j];
                k += 1;
            }
        }
        isTermi[n] = 5;
        isTermi[n+m-1] = 3;
        n += m;
    }
    top->nBP = k;
    // Note!!! top->nBP isn't nBP, the former excluded the lonepair(if minStembp=2)
    

    //-----------------------------------------------------------------------------
    //------------- Step 3. Determine the stacking and angle list -----------------
    //-----------------------------------------------------------------------------
    int a, b, c, d;
    int idx, *anglist = NULL, *stklist = NULL;
    int isStk, nangle, maxNstk, nSTK, iBPstep;
    maxNstk = 4*4*natom;
    anglist = top->anglist;
    stklist = top->stklist;
    for (i = 0; i < maxNstk; i++){
        stklist[i] = 0;
    }
    for (i = 0; i < 2*natom; i++){
        anglist[i] = 0;
    }
    top->nangle = 0;
    nangle = 0;
    nSTK = 0;
    isStk = 0;
    top->nSTk = 0;
    for (i = 0; i < nBP; i++){
        if (isTermi[i]==3){
            // this code deals perfectly with lonepairs, which are always marked as 3'-stem
            isStk = 0;
        }
        if (isStk==1){
            anglist[nangle] = BPlist5[i] * 3;
            anglist[nangle+1] = BPlist3[i] * 3;
            nangle += 2;
            //
            iBPstep = BpStep(BPlist5[i], BPlist5[i+1], BPlist3[i], BPlist3[i+1], top->iac);
            idx = 4*nSTK;
            a = BPlist3[i] * 3;
            b = BPlist5[i] * 3;
            c = BPlist3[i+1] * 3;
            d = BPlist5[i+1] * 3;
            //printf("a=%d, b=%d, c=%d, d=%d\n", a, b, c, d);
            stklist[idx] = a;     // atom-i
            stklist[idx+1] = b;     // atom-j (center)
            stklist[idx+2] = d;     // stom-k
            stklist[idx+3] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+4] = b;     // atom-i
            stklist[idx+5] = a;     // atom-j (center)
            stklist[idx+6] = c;     // stom-k
            stklist[idx+7] = iBPstep;
            stklist[idx+8] = b;     // atom-i
            stklist[idx+9] = d;     // atom-j (center)
            stklist[idx+10] = c;     // stom-k
            stklist[idx+11] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+12] = a;     // atom-i
            stklist[idx+13] = c;     // atom-j (center)
            stklist[idx+14] = d;     // stom-k
            stklist[idx+15] = iBPstep;
            //
            nSTK += 4;
        }
        if (isTermi[i]==5){
            isStk = 1;
            iBPstep = BpStep(BPlist5[i], BPlist5[i+1], BPlist3[i], BPlist3[i+1], top->iac);
            idx = 4*nSTK;
            a = BPlist3[i] * 3;
            b = BPlist5[i] * 3;
            c = BPlist3[i+1] * 3;
            d = BPlist5[i+1] * 3;
            //printf("a=%d, b=%d, c=%d, d=%d\n", a, b, c, d);
            stklist[idx] = a;     // atom-i
            stklist[idx+1] = b;     // atom-j (center)
            stklist[idx+2] = d;     // stom-k
            stklist[idx+3] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+4] = b;     // atom-i
            stklist[idx+5] = a;     // atom-j (center)
            stklist[idx+6] = c;     // stom-k
            stklist[idx+7] = iBPstep;
            stklist[idx+8] = b;     // atom-i
            stklist[idx+9] = d;     // atom-j (center)
            stklist[idx+10] = c;     // stom-k
            stklist[idx+11] = iBPstep;                         // for now, we don't distinguish bp-steps
            stklist[idx+12] = a;     // atom-i
            stklist[idx+13] = c;     // atom-j (center)
            stklist[idx+14] = d;     // stom-k
            stklist[idx+15] = iBPstep;
            //
            nSTK += 4;
        }
    }
    top->nangle = nangle;
    top->nSTk = nSTK;
tim->t_BP_Angs += punchClock(tim->tt);


    //------------- Step 4. Apply harmonic restraint to basepairs ------------------
    double bphm_eq, bphm_k, bphm_p;
    int ntypes = top->pointers[1];
    int idxa, idxb, nHB;
    bphm_eq = top->bphm_eq;
    bphm_k = top->bphm_k;
    bphm_p = top->bphm_p;
    n = 0;
    //printf("Restraint on: ");
    for (i = 0; i < top->nBP; i++){
        // get atom index and force field parm index
        //printf("%d-%d ", BPlist5[i], BPlist3[i]);
        ia = aBPlist5[i]*3;
        ib = aBPlist3[i]*3;
        // get the number of HB for different BP
        idxa = top->iac[aBPlist5[i]]-1;
        idxb = top->iac[aBPlist3[i]]-1;
        iparm = top->ico[ntypes*idxa + idxb] -1;
        nHB = top->nHB[iparm];
        // Compute displacement
        dx = sdt->crd[ia] - sdt->crd[ib];
        dy = sdt->crd[ia+1] - sdt->crd[ib+1];
        //
        r2 =  dx*dx + dy*dy;
        r = sqrt(r2);
        dl = r - bphm_eq;
        //printf("i=%d, ia=%d, ib=%d, dl=%lf\n", i, ia, ib, dl);

        // Accumulate the bond energy and first derivative                       // Notice: it's gradient, not force yet
        *eng += nHB * (bphm_k * dl * dl + bphm_p);
        de = 2.0 * nHB * bphm_k * dl / r;
        deriv[ia] += de*dx;
        deriv[ia+1] += de*dy;
        deriv[ib] -= de*dx;
        deriv[ib+1] -= de*dy;
    }
    //printf("\n");
tim->t_BP_others += punchClock(tim->tt);
    // ------------------------------------------------------------
}




//-----------------------------------------------------------------------------
// eStacking: constrain the 4 rectangular-angles inside the bp-step.
// 
//-----------------------------------------------------------------------------
void eStacking(prmtop *top, state *sdt, double* E_stk){
    int i, j, nangl;
    int ia, ib, ic, iparm, n;
    double eq, k;
    double ba[3], ac[3], bc[3]; // replacement of atom a/b/c at x/y/z axis
    double mgba, mgbc, invbabc;
    double sqba, sqbc, mbabc;
    double costheta, theta, dtheta;
    double de, dea, dec;
    double *eng = NULL;
    double *deriv = NULL; // first derivative
    
    //printf("data has entered eangl()\n");
    //
    eng = E_stk;
    deriv = sdt->acl;

    // get the number of angles (including and without hydrogen)
    int *stklist = NULL, nSTk;
    double stk_p = top->stk_p;
    stklist = top->stklist;
    nSTk = top->nSTk;
    
    //test
    /*
    printf("There are %d Stackings:\n", nangle/2);
    for (i = 0; i < nangle;){
        printf("%d: %d-%d\n", i/2, anglist[i], anglist[i+1]);
        i += 2;
    }
    */

    n = 0;
    eq = top->stk_eq;
    for (i = 0; i < nSTk; i++){
        n = i*4;
        // get atom index and force field parm index
        ia = stklist[n];
        ib = stklist[n+1];
        ic = stklist[n+2];
        iparm = stklist[n+3];

        // get force field parameters;
        k = top->stk_k[iparm];
        // Compute displacement for axis-x/y/z
        for (j = 0; j < 2; j++){
            ba[j] = sdt->crd[j+ia] - sdt->crd[j+ib];
            ac[j] = sdt->crd[j+ic] - sdt->crd[j+ia];
            bc[j] = sdt->crd[j+ic] - sdt->crd[j+ib];
        }

        // On to the angle force computation
        mgba = ba[0]*ba[0] + ba[1]*ba[1];
        mgbc = bc[0]*bc[0] + bc[1]*bc[1];
        invbabc = 1.0/sqrt(mgba*mgbc);
        costheta = (ba[0]*bc[0] + ba[1]*bc[1]) * invbabc;
        //printf("costheta=%lf\t", costheta);
        if (abs(costheta)>0.99){
            continue;
        }
        //costheta = (costheta<-1.0) ? -1.0 : (costheta>1.0) ? 1.0 : costheta;
        theta = acos(costheta);
        ///printf("theta=%lf\n", theta);
        dtheta = theta - eq;
        //printf("eq=%lf, k =%lf, theta =%lf\n", eq, k, theta);

        // Accumulate the angle energy and first derivative
        *eng += k*(dtheta*dtheta + stk_p);
        if (abs(dtheta) < 0.001){
            // no force
            continue;
        }
        de = -2.0*k*dtheta / sqrt(1.0 - costheta*costheta);
        sqba = de/mgba;
        sqbc = de/mgbc;
        mbabc = de * invbabc;

        // loop for axis-x/y/z update of first derivatives
        for (j = 0; j < 2; j++){
            dea = bc[j]*mbabc - costheta*ba[j]*sqba;
            dec = ba[j]*mbabc - costheta*bc[j]*sqbc;
            //printf("dea=%lf, dec=%lf\n", dea, dec);
            deriv[j+ia] += dea;
            deriv[j+ib] -= (dea + dec);
            deriv[j+ic] += dec;
        }
    }
    //printf("Force test after eangl: deriv[0]=%lf\n", deriv[0]);
    //printf("data entered eangl() successfully returned \n");
}
//-----------------------------------------------------------------------------
// eangl: Geometrical restraint of 3 continuous beads in the stem
// 
//-----------------------------------------------------------------------------
void eangl(prmtop *top, state *sdt, double* E_angl){
    int i, j, nangl;
    int ia, ib, ic, iparm, n;
    double eq, k, k_pen;
    double ba[3], ac[3], bc[3]; // replacement of atom a/b/c at x/y/z axis
    double mgba, mgbc, invbabc;
    double sqba, sqbc, mbabc;
    double costheta, theta, dtheta;
    double de, dea, dec;
    double *eng = NULL;
    double *deriv = NULL; // first derivative


    //printf("data has entered eangl()\n");
    //
    eng = E_angl;
    deriv = sdt->acl;

    // get the number of angles (including and without hydrogen)
    int *anglist = NULL, nangle;
    double angle_p = top->angle_p;;
    anglist = top->anglist;
    nangle = top->nangle;    
    
    //test
    /*
    printf("There are %d Stackings:\n", nSTk/2);
    for (i = 0; i < nSTk;){
        printf("%d: %d-%d\n", i/2, stklist[i], stklist[i+1]);
        i += 2;
    }
    */
    // get force field parameters; (uniform stacking constant)
    eq = top->angle_eq;
    k = top->angle_k;
    for (i = 0; i < nangle; i++){
        n = anglist[i];
        // get atom index and force field parm index
        ia = n - 3;
        ib = n;
        ic = n + 3;
        ///printf("ia=%d, ib =%d, ic =%d\n", ia, ib, ic);

        // Compute displacement for axis-x/y/z
        for (j = 0; j < 2; j++){
            ba[j] = sdt->crd[j+ia] - sdt->crd[j+ib];
            ac[j] = sdt->crd[j+ic] - sdt->crd[j+ia];
            bc[j] = sdt->crd[j+ic] - sdt->crd[j+ib];
        }

        // On to the angle force computation
        mgba = ba[0]*ba[0] + ba[1]*ba[1];// + ba[2]*ba[2];
        mgbc = bc[0]*bc[0] + bc[1]*bc[1];// + bc[2]*bc[2];
        invbabc = 1.0/sqrt(mgba*mgbc);  
        costheta = (ba[0]*bc[0] + ba[1]*bc[1]) * invbabc;
        ///printf("costheta=%lf\t", costheta);
        costheta = (costheta<-1.0) ? -1.0 : (costheta>1.0) ? 1.0 : costheta;
        theta = acos(costheta);
        //printf("theta=%lf\n", theta);
        dtheta = theta - eq;
        ////printf("eq=%lf, k =%lf, theta =%lf\n", eq, k, theta);
        // Accumulate the angle energy and first derivative
        *eng += k * dtheta * dtheta + angle_p;
        if (abs(dtheta) < 0.001){/// if theta=pi, the derivative doesn't work!!
            // no force
            continue;
        }
        de = -2.0*k*dtheta / sqrt(1.0 - costheta*costheta);
        sqba = de/mgba;
        sqbc = de/mgbc;
        mbabc = de * invbabc;

        // loop for axis-x/y/z update of first derivatives
        for (j = 0; j < 2; j++){
            dea = bc[j]*mbabc - costheta*ba[j]*sqba;
            dec = ba[j]*mbabc - costheta*bc[j]*sqbc;
            deriv[j+ia] += dea;
            deriv[j+ib] -= (dea + dec);
            deriv[j+ic] += dec;
        }
    }
    //printf("Force test after eangl: deriv[0]=%lf\n", deriv[0]);
    //printf("data entered eangl() successfully returned \n");
}

//-----------------------------------------------------------------------------
// eRestraint: harmonic restriant to designated base-pairs.
// 
//-----------------------------------------------------------------------------
void eRestraint(prmtop *top, state *sdt, double* E_restraint){
    int i, j;
    int ia, ib, iparm, n;
    double eq, k, cutoff;
    double dx, dy, de;
    double r2, r, dl;
    double *eng = NULL;
    double *deriv = NULL;
    int *ntrBplist=NULL, ntr_nBP;

    //
    eng = E_restraint;
    deriv = sdt->acl;

    //
    ntrBplist = top->ntrBplist;
    ntr_nBP = top->ntr_nBP;
    //printf("ntr_nBP=%d\n", ntr_nBP);

    // get force field parameters;
    eq = top->restraint_eq;
    k = top->restraint_wt;
    cutoff = top->rcut * top->rcut;
    //printf("eq=%lf, k =%lf\n", eq, k);

    n = 0;
    for (i = 0; i < ntr_nBP; i++){
        // get atom index and force field parm index
        ia = ntrBplist[n] * 3;
        ib = ntrBplist[n+1] * 3;
        //printf("i=%d, ia=%d, ib =%d\n", i, ia, ib);
        n += 2;
        // Compute displacement
        dx = sdt->crd[ia] - sdt->crd[ib];
        dy = sdt->crd[ia+1] - sdt->crd[ib+1];
        //dz = sdt->crd[ia+2] - sdt->crd[ib+2];
        r2 =  dx*dx + dy*dy;// +dz*dz;
        if (r2>cutoff){
            r = sqrt(r2);
            dl = r - eq;
            //printf("ia=%d, ib=%d, eq=%lf, k =%lf, r =%lf\n", ia, ib, eq, k, r);
            // Accumulate the bond energy and first derivative
            *eng += k*dl*dl;
            de = 2.0*k*dl/r;
            deriv[ia] += de*dx;
            deriv[ia+1] += de*dy;
            deriv[ib] -= de*dx;
            deriv[ib+1] -= de*dy;
        }
    }
}
//-----------------------------------------------------------------------------
// eTargetedMD: the harmonic restraints are applied only after the base-pairs
//              have formed spontaneously.
//-----------------------------------------------------------------------------
void eTargetedMD(prmtop *top, state *sdt, double* E_restraint){
    int i, j, n;
    int ia, ib;
    int ia_rsn, ib_rsn;
    int *BPlist5=NULL, *BPlist3=NULL, nBP;
    int *ntrBplist=NULL, ntr_nBP;
    int *tgtBplist=NULL, *notTargeted=NULL, n_tgtBP;
    int nBP_Left;
    double eq, k, cutoff, r2, r, dl;
    double dx, dy, de;
    double *eng = NULL;
    double *deriv = NULL;
    eng = E_restraint;
    deriv = sdt->acl;


    // get the Basepairs of this conformation (most of them might not be Native)
    nBP = top->nBP;
    BPlist5 = top->BPlist5;     // stored the index of 5'-BP
    BPlist3 = top->BPlist3;     // stored the index of 3'-BP

    // get the designated Basepairs (not necessarily all targeted already)
    ntrBplist = top->ntrBplist;
    ntr_nBP = top->ntr_nBP;

    // get the targeted Basepairs (the BP successfully targeted)
    tgtBplist = top->tgtBplist;
    notTargeted = top->notTargeted;
    n_tgtBP = top->n_tgtBP;

    // check whether all-BP are targeted
    nBP_Left = ntr_nBP - n_tgtBP;
    if (nBP_Left){  // not fully targeted
        // find if any BP are native
        for (i = 0; i < nBP; i++){
            ia = BPlist5[i];
            ib = BPlist3[i];
            for (j = 0; j < ntr_nBP; j++){
                // if this BP[j] is already targeted
                if (notTargeted[j]){
                    ia_rsn = ntrBplist[j*2];
                    ib_rsn = ntrBplist[j*2+1];
                    if (ia==ia_rsn && ib==ib_rsn){
                        tgtBplist[n_tgtBP*2] = ia;
                        tgtBplist[n_tgtBP*2+1] = ib;
                        n_tgtBP += 1;
                        notTargeted[j] = 0;
                    }
                }
            }
        }
    }
    // must store updated n_tgtBP !!!!
    top->n_tgtBP = n_tgtBP;

    // apply the harmonic restraint to the targeted Basepairs
    if (n_tgtBP){
        // get force field parameters;
        eq = top->restraint_eq;
        k = top->restraint_wt;
        cutoff = top->rcut * top->rcut;
        /////printf("eq=%lf, k =%lf， cutoff=%lf\n", eq, k, cutoff);
        n = 0;
        for (i = 0; i < n_tgtBP; i++){
            // get atom index and force field parm index
            ia = tgtBplist[n] * 3;
            ib = tgtBplist[n+1] * 3;
            ////printf("i=%d, ia=%d, ib =%d\n", i, ia, ib);
            n += 2;
            // Compute displacement
            dx = sdt->crd[ia] - sdt->crd[ib];
            dy = sdt->crd[ia+1] - sdt->crd[ib+1];
            //dz = sdt->crd[ia+2] - sdt->crd[ib+2];
            r2 =  dx*dx + dy*dy;// +dz*dz;
            ////printf("r2=%lf\n", r2);
            if (r2>cutoff){
                r = sqrt(r2);
                dl = r - eq;
                ////printf("ia=%d, ib=%d, eq=%lf, k =%lf, r =%lf\n", ia, ib, eq, k, r);
                // Accumulate the bond energy and first derivative
                *eng += k*dl*dl;
                de = 2.0*k*dl/r;
                deriv[ia] += de*dx;
                deriv[ia+1] += de*dy;
                deriv[ib] -= de*dx;
                deriv[ib+1] -= de*dy;
            }
        }
    }
}
//-----------------------------------------------------------------------------
// eSMDcf: Steered MD: constant force stretching
//         (only one pair of residue is allow for now)
//-----------------------------------------------------------------------------
void eSMDcf(prmtop *top, state *sdt, double* E_restraint, double* W_smd){
    int i, j;
    int ia, ib, iparm, n;
    double eq, k;
    double dx, dy, de;
    double r2, r, dl;
    double *eng = NULL, *work = NULL;
    double *deriv = NULL;
    double stretch_force;
    //
    eng = E_restraint;
    work = W_smd;
    deriv = sdt->acl;
    stretch_force = top->stretch_force;

    // get force field parameters;
    // eq = top->restraint_eq;
    k = top->restraint_wt;  // the
    // Work = Force * (X_1 - X_0)
    // X_0 is the old "sdt->spring_len"
    *work -= stretch_force * sdt->spring_len;
    //printf("eq=%lf, k =%lf\n", eq, k);

    // get atom index and force field parm index (default: stretch at the termini)
    ia = 0 * 3;
    ib = (sdt->natom-1) * 3;
    //printf("i=%d, ia=%d, ib =%d\n", i, ia, ib);
    // Compute displacement
    dx = sdt->crd[ia] - sdt->crd[ib];
    dy = sdt->crd[ia+1] - sdt->crd[ib+1];
    //dz = sdt->crd[ia+2] - sdt->crd[ib+2];
    r2 =  dx*dx + dy*dy;// +dz*dz;
    r = sqrt(r2);
    // CF-SMD: change the equil value so that the force is always the designated one
    //  Force = k_spring * (X_t - X_0)
    // X_t = X_0 + Force / k_spring (assume the spring is always slightly longer)
    eq = r + 0.5*stretch_force/k;
    sdt->dist_end = r;
    sdt->spring_len = eq;
    //
    dl = r - eq;
    //printf("ia=%d, ib=%d, eq=%lf, k =%lf, r =%lf\n", ia, ib, eq, k, r);
    // Accumulate the bond energy and first derivative
    *eng += k*dl*dl;
    // Work = Force * (X_1 - X_0)
    // X_1 is the newly updated "sdt->spring_len"
    *work += stretch_force * sdt->spring_len;
    de = 2.0*k*dl/r;
    deriv[ia] += de*dx;
    deriv[ia+1] += de*dy;
    deriv[ib] -= de*dx;
    deriv[ib+1] -= de*dy;
}
//-----------------------------------------------------------------------------
// eSMDcv: Steered MD: constant velocity stretching
//         (only one pair of residue is allow for now)
//-----------------------------------------------------------------------------
void eSMDcv(prmtop *top, state *sdt, double* E_restraint, double* W_smd){
    int i, j;
    int ia, ib, iparm, n;
    double eq, k;
    double dx, dy, de;
    double r2, r, dl;
    double *eng = NULL, *work = NULL;
    double *deriv = NULL;
    //
    eng = E_restraint;
    work = W_smd;
    deriv = sdt->acl;

    // Work = Force * (X_1 - X_0)
    // X_0 is the old "sdt->spring_len"
    *work -= top->stretch_force * sdt->spring_len;
    // get force field parameters;
    // eq = top->restraint_eq;
    k = top->stretch_spring;  // the
    //printf("eq=%lf, k =%lf\n", eq, k);

    // get atom index and force field parm index (default: stretch at the termini)
    ia = 0 * 3;
    ib = (sdt->natom-1) * 3;
    //printf("i=%d, ia=%d, ib =%d\n", i, ia, ib);
    // Compute displacement
    dx = sdt->crd[ia] - sdt->crd[ib];
    dy = sdt->crd[ia+1] - sdt->crd[ib+1];
    //dz = sdt->crd[ia+2] - sdt->crd[ib+2];
    r2 =  dx*dx + dy*dy;// +dz*dz;
    r = sqrt(r2);
    // CF-SMD: change the equil value so that the force is always the designated one
    //  Force = k_spring * (X_t - X_0)
    // X_t = X_0 + Force / k_spring (assume the spring is always slightly longer)
    eq = sdt->spring_len;
    sdt->dist_end = r;
    //
    dl = r - eq;
    // add a 'minus' to the force because we want the force on the pseudoatom
    top->stretch_force = - 2 * k * dl;
    // Work = Force * (X_1 - X_0)
    // X_1 is the newly updated "sdt->spring_len"
    *work += top->stretch_force * sdt->spring_len;
    //printf("ia=%d, ib=%d, eq=%lf, k =%lf, r =%lf\n", ia, ib, eq, k, r);
    // Accumulate the bond energy and first derivative
    *eng += k*dl*dl * r;
    de = 2.0*k*dl/r;
    deriv[ia] += de*dx;
    deriv[ia+1] += de*dy;
    deriv[ib] -= de*dx;
    deriv[ib+1] -= de*dy;
}
//-----------------------------------------------------------------------------
// gradient: computes forces and their first derivatives of all kind in use.
// 
//-----------------------------------------------------------------------------
void gradient(energy *eng, prmtop *top, state *sdt, timer *tim, control *ctl, FILE *fp, int nstep){
    int i, natom;
    natom = sdt->natom;

    // zero out each of the potential energy components
    eng->E_bond = 0.0;          // Energy due to bond-stretching terms
    eng->E_angl = 0.0;
    eng->E_vdw = 0.0;           // Energy due to van der Waarl terms
    eng->E_elec = 0.0;           // Energy due to electrostatic terms
    eng->E_basepair = 0.0;      // Energy due to base-pairing terms
    eng->E_stacking = 0.0;      // Energy due to base-stacking terms
    eng->E_restraint = 0.0;     // Energy due to restraint terms
    eng->Epot = 0.0;            // Total potential energy

    // zero out first derivatives (or 'forces')
    for (i = 0; i < 3*natom; i++){
        sdt->acl[i] = 0.0;
    }
    //-----------------------------------------------------------------------------
    //------------------------ B O N D  -------------------------------------------
    //-----------------------------------------------------------------------------
    ////printf("!! check ebond\n");
    ebond(top, sdt, &eng->E_bond);
    //fprintf(fp, "!! check ebond\n");
    //printf("!!!E_bond=%lf\n", eng->E_bond);
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);

tim->t_Bond += punchClock(tim->tt);


    //-----------------------------------------------------------------------------
    //-------------------- All non-bonded -----------------------------------------
    //-------- Van der Waals + Electrostatics + HB-cmap  --------------------------
    //-----------------------------------------------------------------------------
    eNonBonded(top, sdt, &eng->E_vdw, &eng->E_elec, ctl->dpr, fp, nstep);
    //fprintf(fp, "!! check eNonBonded %d\n", 2);
    //printf("!!!E_vdw=%lf\n", eng->E_vdw);
    //printf("!!!E_elec=%lf\n", eng->E_elec);
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);

    //-----------------------------------------------------------------------------
    //------------------------ Van der Waals --------------------------------------
    //-----------------------------------------------------------------------------
    //evdw(top, sdt, &eng->E_vdw);
    ////printf("!! check evdw %d\n", 2);
tim->t_VDW += punchClock(tim->tt);
    ////printf("!!!E_vdw=%lf\t", eng->E_vdw);
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);

    //-----------------------------------------------------------------------------
    //------------------------ Electrostatic --------------------------------------
    //-----------------------------------------------------------------------------
    //eele(top, sdt, &eng->E_elec);
    ////printf("!! check eele\n");
//tim->t_ELEC += punchClock(tim->tt);
    ////printf("!!!E_elec=%lf\t", eng->E_elec);
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);

    //-----------------------------------------------------------------------------
    //------------------------ Base pairing ---------------------------------------
    //-----------------------------------------------------------------------------
    if (ctl->dpr == 1){
        eBasePair_noDP(top, sdt, &eng->E_basepair, tim);
    }
    // dpr=2 or dpr=0 (non-confined MD)
    else{
        eBasePair(top, sdt, &eng->E_basepair, tim, ctl->dpr, fp, nstep);
    }
    
    //fprintf(fp, "!! check eBasePair\n");
    //printf("!!!E_basepair=%lf\n", eng->E_basepair);
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);

    //-----------------------------------------------------------------------------
    //------------------------ A N G L E  -----------------------------------------
    //-----------------------------------------------------------------------------
    eangl(top, sdt, &eng->E_angl); 
    //printf("!!!E_angl=%lf\n", eng->E_angl);
    //fprintf(fp, "!! check eangl\n");
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);
tim->t_Angle += punchClock(tim->tt);
    

    //-----------------------------------------------------------------------------
    //------------------------ S T A C K I N G  -----------------------------------
    //-----------------------------------------------------------------------------
    eStacking(top, sdt, &eng->E_stacking); 
    //printf("!!!E_stacking is %lf\n", eng->E_stacking);
    //fprintf(fp, "!! check eStacking\n");
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);
tim->t_Stacking += punchClock(tim->tt);
    ////printf("\n");


    //-----------------------------------------------------------------------------
    //------------------------ R E S T R A I N T  ---------------------------------
    //-----------------------------------------------------------------------------
    if (ctl->ntr == 1){
        eRestraint(top, sdt, &eng->E_restraint);
        tim->t_Restraint += punchClock(tim->tt);
        //printf("!!!eRestraint is %lf\n", eng->E_restraint);
    }
    // targeted MD
    else if (ctl->ntr == 2) {
        eTargetedMD(top, sdt, &eng->E_restraint);
        tim->t_Restraint += punchClock(tim->tt);
        //printf("!!!eTargetedMD is %lf\n", eng->E_restraint);
    }
    ////printf("!! check eRestraint\n");
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);


    //-----------------------------------------------------------------------------
    //------------------------ Steered MD  ----------------------------------------
    //-----------------------------------------------------------------------------
    // Steer MD (constant force)
    if (ctl->ismd == 1) {
        eSMDcf(top, sdt, &eng->E_restraint, &eng->W_smd);
        tim->t_Restraint += punchClock(tim->tt);
        //printf("!!!eTargetedMD is %lf\n", eng->E_restraint);
    }
    // Steer MD (constant velocity)
    else if (ctl->ismd == 2) {
        eSMDcv(top, sdt, &eng->E_restraint, &eng->W_smd);
        tim->t_Restraint += punchClock(tim->tt);
        //printf("!!!eTargetedMD is %lf\n", eng->E_restraint);
    }
    ////printf("!! check eSMD\n");
    ////printf("atom-%d has deriv(%d): %.8lf\n", 0, 0, sdt->acl[0]);

    // sum up to get the total energy and first derivatives
    eng->Epot = eng->E_vdw + eng->E_bond + eng->E_basepair + eng->E_angl + eng->E_stacking + eng->E_elec;
}
