#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


//-----------------------------------------------------------------------------
// perfc: this erfc() function was taken from sander source code _erfcfun.f.
//        The algorithm is based on the observation that, for the portion of
//        the domain where erfc() changes rapidly, it appears very much like a
//        product of a Gaussian and a rational function.  Computing these, and
//        optimizing the rational function with seventh-order polynomials,
//        produces the desired result to machine precision.
//
// Arguments:
//   x:  the argument of the complimentary error function
//-----------------------------------------------------------------------------
double perfc(double x)
{
  double absx, c, p, q, nonexperfc, tserf, tserfc;

  absx = fabs(x);
  if (x > 26.0) {
    tserfc = 0.0;
  }
  else if (x < -5.5) {
    tserfc = 2.0;
  }
  else if (absx <= 0.5) {
    c = x*x;
    p = (((-0.356098437018154e-1)*c+6.99638348861914)*c + 21.9792616182942)*c +
      242.667955230532;
    q=((c+15.0827976304078)*c+91.1649054045149)*c + 215.058875869861;
    tserf = x*p/q;
    tserfc = 1.0-tserf;
  }
  else if (absx < 4.0) {
    c = absx;
    p = (((((((-0.136864857382717e-6)*c + 0.564195517478974)*c +
             7.21175825088309)*c + 43.1622272220567)*c +
           152.989285046940)*c + 339.320816734344)*c +
         451.918953711873)*c + 300.459261020162;
    q = ((((((c+12.7827273196294)*c + 77.0001529352295)*c +
            277.585444743988)*c + 638.980264465631)*c +
          931.354094850610)*c + 790.950925327898)*c + 300.459260956983;
    nonexperfc = (x > 0.0) ? p/q : 2.0*exp(x*x) - p/q;
    tserfc = exp(-absx*absx)*nonexperfc;
    if (x < 0.0) {
      tserfc = 2.0 - tserfc;
    }
  }
  else {
    c = 1.0/(x*x);
    p=(((0.0223192459734185*c + 0.278661308609648)*c +
        0.226956593539687)*c + 0.0494730910623251)*c + 0.00299610707703542;
    q=(((c + 1.98733201817135)*c + 1.05167510706793)*c + 0.191308926107830)*c +
      0.0106209230528468;
    c = (-c*p/q + 0.564189583547756)/absx;
    nonexperfc = (x > 0.0) ? c : 2.0*exp(x*x) - c;
    tserfc = exp(-absx*absx)*nonexperfc;
    if (x < 0.0) {
      tserfc = 2.0 - tserfc;
    }
  }

  return tserfc;
}

//-----------------------------------------------------------------------------
// PerfcFast: this is a faster way to compute the complimentary error function.
//            It assumes that x is positive, which is usually the case in the
//            sorts of calculations done for molecular dynamics simulations.
//                                                                      
// Arguments:                                                           
//   x:  the argument of the complimentary error function               
//-----------------------------------------------------------------------------
double PerfcFast(double x)
{
  double c, p, q, tserf, tserfc, invx;

  // Most likely, x will be between 0.5 and 4.0
  if (x > 0.5 && x < 4.0) {
    p = (((((((-0.136864857382717e-6)*x + 0.564195517478974)*x +
             7.21175825088309)*x + 43.1622272220567)*x +
           152.989285046940)*x + 339.320816734344)*x +
         451.918953711873)*x + 300.459261020162;
    q = ((((((x+12.7827273196294)*x + 77.0001529352295)*x +
            277.585444743988)*x + 638.980264465631)*x +
          931.354094850610)*x + 790.950925327898)*x + 300.459260956983;
    tserfc = exp(-x*x)*p/q;
    return tserfc;
  }

  // In a minority of cases, x will be less than 0.5
  if (x <= 0.5) {
    c = x*x;
    p = (((-0.356098437018154e-1)*c+6.99638348861914)*c + 21.9792616182942)*c +
      242.667955230532;
    q=((c+15.0827976304078)*c+91.1649054045149)*c + 215.058875869861;
    tserf = x*p/q;
    tserfc = 1.0-tserf;
    return tserfc;
  }

  // Otherwise, we have one more condition to evaluate about x
  if (x >= 4.0 && x < 5.65) {
    invx = 1.0/x;
    c = invx*invx;
    p=(((0.0223192459734185*c + 0.278661308609648)*c +
        0.226956593539687)*c + 0.0494730910623251)*c + 0.00299610707703542;
    q=(((c + 1.98733201817135)*c + 1.05167510706793)*c + 0.191308926107830)*c +
      0.0106209230528468;
    c = (-c*p/q + 0.564189583547756)*invx;
    tserfc = exp(-x*x)*c;
    return tserfc;
  }
  
  // If x is bigger than 5.65, the answer is easy
  return 0.0;
}

//-------------------------------------------------------------------------
// CrossP: function for finding the cross-product cr of vectors p and q.  
//         Note that vectors p and q are assumed to be three-dimensional 
//         and only the first three components of these vectors will be 
//         considered.
//-------------------------------------------------------------------------
void CrossP(double* p, double* q, double* cr){
    cr[0] = p[1]*q[2]-p[2]*q[1];
    cr[1] = p[2]*q[0]-p[0]*q[2];
    cr[2] = p[0]*q[1]-p[1]*q[0];
}

//-----------------------------------------------------------------------------
// UnitVector3: normalize a vector V of length 3.                       
//-----------------------------------------------------------------------------
void UnitVector3(double* V){
    double mg;

    mg = 1.0/sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
    V[0] *= mg;
    V[1] *= mg;
    V[2] *= mg;
}

//-----------------------------------------------------------------------------
// Dot3: Dot product of two double-precision real vectors V1 and V2 of length
//       3.  Useful for coordinate or displacement dot products. 
//-----------------------------------------------------------------------------
double Dot3(double* V1, double* V2)
{
  return V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
}