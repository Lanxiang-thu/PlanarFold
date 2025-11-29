#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dsfly.h"

//-----------------------------------------------------------------------------
// RAN2: function for returning a single random number from a uniform   
//       distribution in the range (0, 1), exclusive of the endpoints.  
//                                                                      
// Arguments:                                                           
//   nwin:  pointer to a long unsigned integer, the window into the 3-digit
//          random number state vector                         
//-----------------------------------------------------------------------------
double ran2(long *nwin){ 
  	int j;
  	long k, nwinlocal;
  	static long nwin2=123456789;
  	static long iy=0;
  	static long iv[32];
  	double temp;

  	// Local value for the modifiable input variable
  	nwinlocal = *nwin;

	// Initialize. Be sure to prevent nSeed = 0
  	if (nwinlocal <= 0) {
    	nwinlocal = (-nwinlocal < 1) ? 1 : -nwinlocal;
    	nwin2 = nwinlocal;//nwin2=1
		// Load the shuffle table ( after 8 warm-ups )
    	for (j = 32+7; j >= 0 ; j--) {
      		k = nwinlocal/53668;
      		nwinlocal = 40014*(nwinlocal - k*53668) - k*12211;
      		if (nwinlocal < 0) {
				nwinlocal += 2147483563;
      		}
      		if (j < 32) {
				iv[j] = nwinlocal;
      		}
    	}
    	iy = iv[0];
  	}

	// start here when not initializing.
  	k = nwinlocal/53668;

	// Compute idum=(IA*idum)%IM without overflows by Schrage's method
  	nwinlocal = 40014*(nwinlocal - k*53668) - k*12211; 
  	if (nwinlocal < 0) {
    	nwinlocal += 2147483563;
  	}
  	k = nwin2/52774;
  	nwin2 = 40692*(nwin2 - k*52774) - k*3791; 
  	if (nwin2 < 0) {
    	nwin2 += 2147483399;
  	}
	// Will be in the range 0..NTAB-1.
  	j = iy / (1 + 67108861);// (int) (2147483562/32) = 67108861

	// Output previously stored value and refill the shuffle table
  	iy = iv[j] - nwin2; 
  	iv[j] = nwinlocal;

  	// Copy local value back to modifiable input variable
  	*nwin = nwinlocal;
  	if (iy < 1) iy += 2147483562;// ~ 2^31
  	temp = iy * 1.0 / 2147483563;
  	if (temp > 0.99999988) {// 1.0-1.2e-7
    	return 0.99999988;
  	}
    /*else if (temp < 0.00000001) {// 1.0-1.2e-7
    	return 0.00000001;
  	}
    */
  	else {
    	return temp;
  	}
}
//-----------------------------------------------------------------------------
// Gaussian_random: the vel obeys gausian distribution
//-----------------------------------------------------------------------------
double generate_gaussian_random(double u1, double u2) {
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void verlet(control *ctl, prmtop *top, energy *eng,
    state *sdt, ncfile *ncf, timer *tim, int nstep, FILE *fp){
    int i, j, k, natom, i3;
    double invmass, mass;
    const double dt = ctl->dt;  // sqrt(418.4) or just 418.4??              Only weird mdgx codes will add 418.4 here
    const double hdt = 0.5*dt;

    //printf("data has entered verlet()\n");
    natom = sdt->natom;

    // store the current atom coordinates, then find half-step
    // velocities and full-step coordinates via Verlet recursion
    i3 = 0;
    //printf("Full-step coordinates\n");
    for (i = 0; i < natom; i++){
        //printf("%d, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf\n", i, sdt->crd[i3], sdt->crd[i3+1], sdt->vel[i3], sdt->vel[i3+1], sdt->acl[i3], sdt->acl[i3+1]);
        for (j = 0; j < 2; j++){
            // get half-step velocities
            sdt->vel[i3+j] += sdt->acl[i3+j]*hdt;
            // get the full-step coordinates
            sdt->crd[i3+j] += sdt->vel[i3+j]*dt;
        }
        //printf("%d, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf, %.8lf\n", i, sdt->crd[i3], sdt->crd[i3+1], sdt->vel[i3], sdt->vel[i3+1], sdt->acl[i3], sdt->acl[i3+1]);
        i3 += 3;
    }
tim->t_Verlet += punchClock(tim->tt);

    // update the verlet list (every nVerlet steps)
    if ((nstep-1) % ctl->nVerlet == 0){
        int ia3, ib3, n;
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
    }

    //  get the potential energy and update atomic forces at full step
    //printf("about to calculate the gradient\n");
    gradient(eng, top, sdt, tim, ctl, fp, nstep);
    //printf("the gradient() is okay!\n");

    // use Newton second law to get the next accelerations;
    // find the full-step velocities using the Verlet recursion
    i3 = 0;
    //printf("Full-step velocities\n");
    for (i = 0; i < natom; i++){
        invmass = top->inv_mass[i] * -418.4;           // 418.4 will convert kcal to J
        for (j = 0; j < 2; j++){
            // now we get the real 'acceleration' from former cumputed
            // first derivatives in gradient()
            sdt->acl[i3+j] *= invmass;
            // get the full-step velocities
            sdt->vel[i3+j] += sdt->acl[i3+j]*hdt;
        }
        //printf("%d, %lf, %lf\t", i, sdt->vel[i3], sdt->vel[i3+1]);
        i3 += 3;
    }
tim->t_Verlet += punchClock(tim->tt);

    // Compute kinetic energy and instantaneous temperature after the force computation  
    double velSum;
    double tscale;
    double lambda;
    i3 = 0;
    // kinetic energy: Ek
    eng->Ekin = 0.0;    // zero out kinetic energy
    for (i = 0; i < natom; i++){
        /*
        for (j = 0; j < 2; j++){
            if (isnan(sdt->vel[i3+j])){
                printf("atom-%d has velociy(%d): %lf\n", i, j, sdt->vel[i3+j]);
                exit(0);
            }
        }*/
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
    /*
    printf("temperature=%lf, Ektot=%lf, Eptot=%lf\n", sdt->T, eng->Ekin, eng->Epot);
    if (fabs(sdt->T) > 2000.0){
        printf("temperature %lf larger than 2000.0K \n", sdt->T);
        i3 = 0;
        for (i = 0; i < natom; i++){
            for (j = 0; j < 2; j++){
                if (fabs(sdt->vel[i3+j]) > 410.0){
                    printf("atom-%d has velociy(%d): %lf\n", i, j, sdt->vel[i3+j]);
                    exit(0);
                }
            }
            i3 += 3;
        }
        exit(0);
    }
    */

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
     
        //tscale = (sdt->Ttarget / sdt->T) - 1.0;
        // get the scaling factor 'lambda'
        lambda = sqrt(1.0 + tscale * dt / ctl->tautp);
        // velocity scaling (no virtual sites)
        for (i = 0; i < 3*natom; i++){
            sdt->vel[i] *= lambda;
        }
    }
    /*
    printf("Full-step velocities (after thermostat)\n");
    for (i = 0; i < natom; i++){
        printf("%d, %lf, %lf\t", i, sdt->vel[i*3], sdt->vel[i*3+1]);
    }*/
tim->t_Thermostat += punchClock(tim->tt);

    //-------- Remove momentum (both Translation and Rotation)--------------------------------------------
    double mvx, mvy, mvz;
    double Ltot, dL, rx, ry, J, w;
    double crdCOM[2];// cosT, sinT, dist, vel;
    if (nstep % ctl->nscm == 0){  //i>0 cause we don't wish crd[0] to be stored
        //=======================================================================
        // Random collision of water (***** important for long RNA  ********)
        //=======================================================================
        double dice1, dice2, ratio;    // the random number for movement determination, just like the 'dice'
        double velWat;  // the velocity if the simulated water
        // the collision is not perfectly elastic collision
        double coefR = 0.2; //coefficient of restitution
        if (ctl->iCollide == 1){
            // asign a collision into each residue, mimicking water collision
            for (i = 0; i < natom; i++){
                i3 = 0;
                for (j = 0; j < 2; j++){
                    // dice1 and dice2 are uniformly distributed random number
                    dice1 = ran2(&sdt->rndcon);
                    dice2 = ran2(&sdt->rndcon);
                    // ratio is the Gaussian distributed random number
                    ratio = generate_gaussian_random(dice1, dice2);
                    // velWat = 0.0 + random * standard-deviation
                    // Standard-deviation: from Gaussian distribution
                    //                =  kB * T * NA / mass-WAT / J-to-kcal
                    velWat = ratio * sqrt(1.38 * sdt->T * 6.02 / 18.0 / 4184);
                    //printf("dice1=%lf, dice2=%lf,ratio=%lf, velWat=%lf\n", dice1, dice2, ratio, velWat);
                    // the collision is not perfectly elastic collision
                    sdt->vel[i3+j] += velWat * coefR;
                    //sdt->vel[i3+2] -= mvz;
                }
                i3 += 3;
            }
        }
        //-----------------------------------------------------------------------------
        // remove momentum
        mvx = 0.0;
        mvy = 0.0;
        //mvz = 0.0;
        i3 = 0;
        for (i = 0; i < natom; i++){
            mass = top->mass[i];
            mvx += sdt->vel[i3] * mass;
            mvy += sdt->vel[i3+1] * mass;
            //mvz += sdt->vel[i3+2] * mass;
            i3 += 3;
        }
        mvx /= top->MassTot;
        mvy /= top->MassTot;
        //mvz /= top->MassTot;
        i3 = 0;
        for (i = 0; i < natom; i++){
            sdt->vel[i3] -= mvx;
            sdt->vel[i3+1] -= mvy;
            //sdt->vel[i3+2] -= mvz;
            i3 += 3;
        }
        //-----------------------------------------------------------------------------
        // remove angular momentum
        // locate the Center-Of-Mass
        crdCOM[0] = 0.0;
        crdCOM[1] = 0.0;
        k = 0;
        for (i = 0; i < natom; i++){
            mass = top->mass[i];
            crdCOM[0] += sdt->crd[k] * mass;
            crdCOM[1] += sdt->crd[k+1] * mass;
            k += 3;
        }
        crdCOM[0] /= top->MassTot;
        crdCOM[1] /= top->MassTot;
        // calculate the angular momentum of the system
        Ltot = 0.0;
        J = 0.0;
        k = 0;
        for (i = 0; i < natom; i++){
            mass = top->mass[i];
            rx = sdt->crd[k] - crdCOM[0];
            ry = sdt->crd[k+1] - crdCOM[1];
            Ltot += mass * (rx * sdt->vel[k+1] - ry * sdt->vel[k]);
            J += mass * (rx*rx + ry*ry);
            k += 3;
        }
        // angular velocity = (angular momentum) / (angular inertia)
        w = Ltot / J;
        //printf("Ltot=%lf, J=%lf\t", Ltot, J);
        // remove rotation
        k = 0;
        for (i = 0; i < natom; i++){
            mass = top->mass[i];
            rx = sdt->crd[k] - crdCOM[0];
            ry = sdt->crd[k+1] - crdCOM[1];
            //dist = sqrt(rx*rx + ry*ry);
            //vel = sqrt(sdt->vel[k]*sdt->vel[k] + sdt->vel[k+1]*sdt->vel[k+1]);
            //cosT = fabs((rx*sdt->vel[k] + ry*sdt->vel[k+1]) / (dist * vel));
            //sinT = sqrt(1 - cosT*cosT);
            //printf("cosT=%lf, sinT=%lf\t", cosT, sinT);
            sdt->vel[k] += w * ry;      ///!!!!!
            sdt->vel[k+1] -= w * rx;
            k += 3;
        }
        //printf("\n");
    }
/*
    //--------------------- recalculate kinetic energy ------------------------------------
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
*/
tim->t_Thermostat += punchClock(tim->tt);
}

