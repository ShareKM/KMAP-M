/*
 * kinlib.cpp
 * 
 * This file contains the core kinetic modeling functions used to calculate time activity curves (TAC) and sensitivity functions for different kinetic models. 
 * The Levenberg-Marquardt optimization algorithm is also implemented for solving least squares problems with constraints.
 */

#include "kinlib.h"
#include "mex.h"

//------------------------------------------------------------------------------
// tac_eval
//------------------------------------------------------------------------------
// Evaluates the time activity curve (TAC) for the provided kinetic model parameters.
void tac_eval(double *p, void *param, double *tac)
{ 
   KMODEL_T *par;
   par = (KMODEL_T *) param;
   (*(par->tacfunc))(p, par->dk, par->scant, par->td, par->cp, par->wb, 
                     par->num_frm, par->num_vox, tac);
}

//------------------------------------------------------------------------------
// jac_eval
//------------------------------------------------------------------------------
// Evaluates the Jacobian (sensitivity matrix) for the provided kinetic model parameters.
void jac_eval(double *p, void *param, double *tac, int *psens, double *jac)
{
   KMODEL_T *par;
   par = (KMODEL_T *) param;
   (*(par->jacfunc))(p, par->dk, par->scant, par->td, par->cp, par->wb, 
                     par->num_frm, par->num_vox, tac, psens, jac);
}

//------------------------------------------------------------------------------
// kconv_2t5p_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the two-tissue compartment model.
// This model uses five parameters (2T5P): vb, k1, k2, k3, and k4.
void kconv_2t5p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, k3, k4;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t;
   int      num_time;

   // Memory allocation for intermediate calculations
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);

   for (j=0; j<num_vox; j++) {
      // Parameter transformations
      vb = p[0+j*5];
      k1 = p[1+j*5];
      k2 = p[2+j*5];
      k3 = p[3+j*5];
      k4 = p[4+j*5];
      k234 = k2+k3+k4;
      d = sqrt(k234*k234 - 4*k2*k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;

      // Convolution with exponential functions
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // Compute tissue concentration
      if (d==0) d = 1e9;
      f1 = k1/d * (k4 - a1);
      f2 = k1/d * (a2 - k4);
      b1 = k1/d * k3;
      b2 = -b1;

      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * wb[i];
      }

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
}

//------------------------------------------------------------------------------
// kconv_2t5p_jac
//------------------------------------------------------------------------------
// Calculate the time activity curves and the sensitivity functions 
// (Jacobian) for the two-tissue compartment model.
void kconv_2t5p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, k3, k4;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   int      num_time;
   int      num_par;
  
   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 5; i++) {
      if (psens[i] == 1)
         ++num_par;
   }
  
   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);
   c_t  = (double*) malloc(sizeof(double) * num_time);
   s_t  = (double*) malloc(sizeof(double) * num_time * num_par);
	
   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {
	
      // Transform kinetic parameters into exponential parameters
      vb = p[0 + j * 5];
      k1 = p[1 + j * 5];
      k2 = p[2 + j * 5];
      k3 = p[3 + j * 5];
      k4 = p[4 + j * 5];
      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;
	
      // Calculate the exponential components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // Calculate concentrations
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;	
      for (i = 0; i < num_time; i++) {
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * wb[i];
      }

      // Sensitivity functions (Jacobian) computation
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + wb[i];
         s_t += num_time;
      }
      if (psens[1] == 1 || psens[2] == 1) { 
         f1 = 1 / d * (k4 + k3 - a1);
         f2 = 1 / d * (a2 - k4 - k3);
      }
      if (psens[1] == 1) { // Sensitivity wrt k1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[2] == 1 || psens[3] == 1) { 
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a2);
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[3] == 1 || psens[4] == 1) { 
         f1 = 1 / d * (a1 + a2 - k3 - k4);
         f2 = -f1;
      }
      if (psens[3] == 1) { // Sensitivity wrt k3
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[4] == 1) { // Sensitivity wrt k4
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a2);
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Frame-averaged activity and sensitivity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }
  	
   // Free allocated memory
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------
// kconv_srtm_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the simplified reference tissue model (SRTM).
// This model uses four parameters: vb, R1, k2, and BP.
void kconv_srtm_tac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, R1, k2, BP;
   double   k2a, a;
   double  *c_a, *c_t;
   int      num_time;

   // Memory allocation for intermediate calculations
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);

   for (j = 0; j < num_vox; j++) {
      // Parameter transformations
      vb = p[0 + j * 4];
      R1 = p[1 + j * 4];
      k2 = p[2 + j * 4];
      BP = p[3 + j * 4];

      // Convolution with exponential functions
      k2a = k2 / (1.0 + BP);
      a = k2a + dk;
      kconv_exp(1.0, a, cr0, num_time, td, c_a);

      // Compute tissue concentration
      for (i = 0; i < num_time; i++)
         c_t[i] = (1 - vb) * (R1 * cr0[i] + (k2 - R1 * k2a) * c_a[i]) + vb * wb[i];

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
}

//------------------------------------------------------------------------------
// kconv_srtm_jac
//------------------------------------------------------------------------------
// Calculates the time activity curves and sensitivity functions for the 
// Simplified Reference Tissue Model (SRTM).
void kconv_srtm_jac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, R1, k2, BP;
   double   k2a;
   double   *c_a, *c_b, *c_t, *s_t;
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 4; i++) { 
      if (psens[i] == 1)
         ++num_par;
   }

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a = (double*) malloc(sizeof(double) * num_time);
   c_b = (double*) malloc(sizeof(double) * num_time);
   c_t = (double*) malloc(sizeof(double) * num_time);
   s_t = (double*) malloc(sizeof(double) * num_time * num_par);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Assign parameters
      vb = p[0 + j * 4];
      R1 = p[1 + j * 4];
      k2 = p[2 + j * 4];
      BP = p[3 + j * 4];

      // Compute k2a
      k2a = k2 / (1.0 + BP);
      
      // Calculate the exponential components
      kconv_exp(1.0, k2a + dk, cr0, num_time, td, c_a);

      // Calculate concentration c_t
      for (i = 0; i < num_time; i++) {
         c_b[i] = R1 * cr0[i] + (k2 - R1 * k2a) * c_a[i];
         c_t[i] = (1 - vb) * c_b[i] + vb * wb[i];
      }

      // Average the values over frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);

      // Calculate sensitivity functions for each parameter
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = -c_b[i] + wb[i];
         s_t += num_time;
      }
      if (psens[1] == 1) { // Sensitivity wrt R1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (cr0[i] - k2a * c_a[i]);
         s_t += num_time;
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         for (i = 0; i < num_time; i++)
            c_t[i] = cr0[i] - c_b[i] / (1.0 + BP);
         kconv_exp(1.0, k2a + dk, c_t, num_time, td, c_a);
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * c_a[i];
         s_t += num_time;
      }
      if (psens[3] == 1) { // Sensitivity wrt BP
         kconv_exp(k2a / (1 + BP), k2a + dk, c_b, num_time, td, c_a);
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * c_a[i];
         s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Average the sensitivity functions over frame duration
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }

   // Free allocated memory
   free(c_a);
   free(c_b);
   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------
// kconv_1t3p_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the one-tissue compartment model.
// This model uses three parameters: vb, k1, and k2.
void kconv_1t3p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2;
   double   a;
   double  *c_a, *c_t;
   int      num_time;

   // Memory allocation
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);

   for (j = 0; j < num_vox; j++) {
      // Assign parameters
      vb = p[0 + j * 3];
      k1 = p[1 + j * 3];
      k2 = p[2 + j * 3];

      // Convolution with exponential functions
      a = k2 + dk;
      kconv_exp(k1, a, cp, num_time, td, c_a);

      // Compute tissue concentration
      for (i = 0; i < num_time; i++)
         c_t[i] = (1 - vb) * c_a[i] + vb * wb[i];

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
}

//------------------------------------------------------------------------------
// kconv_1t3p_jac
//------------------------------------------------------------------------------
// Calculates the time activity curves and the sensitivity functions for the 
// one-tissue kinetic model.
void kconv_1t3p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2;
   double   a;
   double  *c_a, *c_t, *s_t;
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 3; i++) {
      if (psens[i] == 1) ++num_par;
   }

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a = (double*) malloc(sizeof(double) * num_time);
   c_t = (double*) malloc(sizeof(double) * num_time);
   s_t = (double*) malloc(sizeof(double) * num_time * num_par);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Assign parameters
      vb = p[0 + j * 3];
      k1 = p[1 + j * 3];
      k2 = p[2 + j * 3];

      // Calculate the exponential component
      a = k2 + dk;
      kconv_exp(k1, a, cp, num_time, td, c_a);

      // Calculate the time activity curve
      for (i = 0; i < num_time; i++)
         c_t[i] = (1 - vb) * c_a[i] + vb * wb[i];

      // Calculate sensitivity functions
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = -c_a[i] + wb[i];
         s_t += num_time;
      }
      if (psens[1] == 1) { // Sensitivity wrt k1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * c_a[i] / k1;
         s_t += num_time;
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         kconv_exp(1, a, c_a, num_time, td, s_t);
         for (i = 0; i < num_time; i++)
            s_t[i] *= -(1 - vb);
         s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Average the sensitivity functions over the frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
   free(s_t);
}



//------------------------------------------------------------------------------
// frame
//------------------------------------------------------------------------------
// Averages the time activity curve (TAC) over the frame duration.
// This function calculates the average activity in the specified time frames.
void frame(double *scant, double td, double *c_t, int num_frm, int num_c, double *ct)
{
   int      m, k;
   double   tmp, dt;
   int      num_time;
   int      j, jstart, jend;

   num_time = (int) (scant[2*num_frm-1] / td);

   for (m = 0; m < num_frm; m++) {
      jstart = (int) (scant[m] / td);
      jend = (int) (scant[m + num_frm] / td);
      dt = (double) (jend - jstart);
      for (k = 0; k < num_c; k++) {
         tmp = 0.0;
         for (j = jstart; j < jend; j++) {
            tmp += c_t[j + k * num_time];
         }
         ct[m + k * num_frm] = tmp / dt;
      }
   }
}

//------------------------------------------------------------------------------
// kconv_exp
//------------------------------------------------------------------------------
// Computes the convolution of an exponential function with the input function.
// c(t) = k1 * exp(-k2 * t) ⊗ u(t) where ⊗ represents the convolution operator.
void kconv_exp(double k1, double k2, double *u, int num_time, double td, double *c)
{
   int      i;
   double   ek2, prev, tmp;

   k1 *= td / 60.0;  // Convert td from seconds to minutes
   k2 *= td / 60.0;
   ek2 = exp(-k2);
   prev = 0;

   if (ek2 == 1.0) {  // Special case when k2 is very small
      for (i = 0; i < num_time; i++) {
         prev += u[i];
         c[i] = k1 * prev;
      }
   } else {
      tmp = (1.0 - ek2) / k2;
      for (i = 0; i < num_time; i++) {
         prev = prev * ek2 + k1 * tmp * u[i];
         c[i] = prev;
      }
   }
}

//------------------------------------------------------------------------------
// kmap_levmar
//------------------------------------------------------------------------------
// The Levenberg-Marquardt algorithm to solve the nonlinear least squares problem
// subject to linear bounds.
// This function is used to estimate the parameters that best fit the model to the data.
void kmap_levmar(double *y, double *w, int num_frm, double *pinit, int num_p, 
                 void *param, void (*func)(double *, void *, double *),
                 void (*jacf)(double *, void *, double *, int *, double *),
                 double *plb, double *pub, int *psens, int maxit, double *ct)
{
   int      i, j, it, n;
   double   mu, v, tau, rho, maxaii;
   int      num_par;
   double   *f, *fnew, *r;
   double   *h;
   double   *pnew, *pest, *plb_sens, *pub_sens;
   double   *st;
   double   etol = 1e-9;
   double   F, Fnew, Lnew;
   int      subit = 100;
   double   tmp;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (j = 0; j < num_p; j++)
      if (psens[j] == 1) ++num_par;

   // Memory allocation
   f = (double*) malloc(sizeof(double) * num_frm);
   fnew = (double*) malloc(sizeof(double) * num_frm);
   r = (double*) malloc(sizeof(double) * num_frm);
   h = (double*) malloc(sizeof(double) * num_par);
   pnew = (double*) malloc(sizeof(double) * num_par);
   pest = (double*) malloc(sizeof(double) * num_par);
   plb_sens = (double*) malloc(sizeof(double) * num_par);
   pub_sens = (double*) malloc(sizeof(double) * num_par);
   st = (double*) malloc(sizeof(double) * num_frm * num_par);

   // Initialize the parameter estimates and their bounds
   getkin(pinit, psens, num_par, pest);
   getkin(plb, psens, num_par, plb_sens);
   getkin(pub, psens, num_par, pub_sens);

   // Initial TAC and sensitivity matrix
   (*jacf)(pinit, param, ct, psens, st);

   // Compute the initial residual and the cost value
   for (i = 0; i < num_frm; i++) {
      f[i] = y[i] - ct[i];
   }
   F = vecnormw(w, f, num_frm) * 0.5;

   // Initialize parameters for the search process
   v = 2.0;
   tau = 1.0e-3;
   maxaii = 0;
   for (j = 0; j < num_par; j++) {
      tmp = vecnorm2(st + num_frm * j, num_frm);
      tmp *= tmp;
      if (maxaii < tmp)
         maxaii = tmp;
   }
   mu = tau * maxaii;

   // Iterative optimization loop
   it = 0; n = 0;
   while (it < maxit) {
      // Estimate new parameters
      for (j = 0; j < num_par; j++)
         pnew[j] = pest[j];

      // Coordinate least squares step
      boundpls_cd(f, w, st, mu, num_frm, num_par, plb_sens, pub_sens, subit, pnew, r);

      // Check for convergence or update mu
      for (i = 0; i < num_par; i++) {
         h[i] = pnew[i] - pest[i];
      }
      if (vecnorm2(h, num_par) <= etol * (vecnorm2(pest, num_par) + etol))
         break;

      // Calculate the ratio rho
      setkin(pnew, num_par, psens, pinit);
      (*func)(pinit, param, ct);
      for (i = 0; i < num_frm; i++)
         fnew[i] = y[i] - ct[i];
      Fnew = vecnormw(w, fnew, num_frm) * 0.5;
      Lnew = vecnormw(w, r, num_frm) * 0.5;
      rho  = (F - Fnew) / (F - Lnew);

      // Accept the step or find a new one
      if (rho > 0) { 
         // Update estimates
         for (j = 0; j < num_par; j++) {
            pest[j] = pnew[j];
         }
         (*jacf)(pinit, param, ct, psens, st);
         for (i = 0; i < num_frm; i++) {
            f[i] = fnew[i];
         }
         F = Fnew;
         ++it;

         // Update mu
         tmp = 2 * rho - 1;
         tmp = 1 - tmp * tmp * tmp;
         if (tmp < 1.0 / 3.0) // must be 1.0/3.0, 1/3 doesn't work
            tmp = 1.0 / 3.0;
         mu *= tmp;
         v = 2.0;
      } else {  
         // Find a new step
         mu *= v;
         v *= 2.0;
      }
      ++n;
      if (n > 1000) {
         printf("stopped: maximum iteration number exceeds 1000!\n");
         setkin(pest, num_par, psens, pinit);
         (*func)(pinit, param, ct);
         break;
      }
   }

   // Free allocated memory
   free(f);
   free(fnew);
   free(r);
   free(h);
   free(pnew);
   free(pest);
   free(plb_sens);
   free(pub_sens);
   free(st);
}

//------------------------------------------------------------------------------
// boundpls_cd
//------------------------------------------------------------------------------
// Coordinate descent algorithm to solve the penalized least squares problem
// subject to box bounds. Minimizes || y - A * (x - x_n) ||^2.
void boundpls_cd(double *y, double *w, double *a, double alpha, int num_y,
                 int num_x, double *xlb, double *xub, int maxit, double *x, 
                 double *r)
{
   int      i, j, ij; 
   double   err0, err, dj, xj;
   double   tmp1, tmp2;
   double   etol = 1.0e-9;
   int      it;

   // Initialize the residual
   for (i = 0; i < num_y; i++){
      r[i] = y[i];
   }
   err0 = vecnorm2(y, num_y);

   // Iterative optimization
   for (it = 0; it < maxit; it++) {
      // Check for convergence
      err = vecnorm2(r, num_y);
      if (err / err0 < etol)
         break;

      // Coordinate-wise optimization
      for (j = 0; j < num_x; j++) {
         tmp1 = tmp2 = 0;
         for (i = 0; i < num_y; i++) {
            ij = i + j * num_y;
            tmp1 += w[i] * a[ij] * r[i];
            tmp2 += w[i] * a[ij] * a[ij];
         }
         tmp2 += alpha * tmp2;
         if (tmp2 == 0)
            dj = 0;
         else
            dj = tmp1 / tmp2;
         xj = x[j] + dj;

         // Apply bounds to the estimate
         if (xj < xlb[j])
            xj = xlb[j];
         else if (xj > xub[j])
            xj = xub[j];

         // Update x and r
         dj = xj - x[j];
         x[j] = xj;
         for (i = 0; i < num_y; i++)
            r[i] -= a[i + j * num_y] * dj;
      }
   }
}

//------------------------------------------------------------------------------
// setkin
//------------------------------------------------------------------------------
// Assigns sensitive parameters from x to x0 based on the sensitivity indicators.
void setkin(double *x, int num_x, int *xsens, double *x0)
{
   int j = 0;
   int k = 0;
   while (j < num_x) {
      if (xsens[k] == 0)
         ++k;
      else {
         x0[k] = x[j];
         ++j;
         ++k;
      }
   }
}

//------------------------------------------------------------------------------
// getkin
//------------------------------------------------------------------------------
// Retrieves sensitive parameters from x0 to x based on the sensitivity indicators.
void getkin(double *x0, int *xsens, int num_x, double *x)
{
   int j = 0;
   int k = 0;
   while (j < num_x) {
      if (xsens[k] == 0)
         ++k;
      else {
         x[j] = x0[k];
         ++j;
         ++k;
      }
   }
}

//------------------------------------------------------------------------------
// vecnorm2
//------------------------------------------------------------------------------
// Computes the L2 norm (Euclidean norm) of vector x.
double vecnorm2(double *x, int num)
{
   int      i;
   double   temp = 0.0;  
  
   for (i = 0; i < num; i++)
      temp += x[i] * x[i];
   return sqrt(temp);
}

//------------------------------------------------------------------------------
// vecnormw
//------------------------------------------------------------------------------
// Computes the weighted norm of vector x.
double vecnormw(double *w, double *x, int num)
{
   int      i;
   double   temp = 0.0;  
  
   for (i = 0; i < num; i++)
      temp += w[i] * x[i] * x[i];
   return temp;
}

//------------------------------------------------------------------------------
// BoundQuadCD
//------------------------------------------------------------------------------
// Coordinate descent algorithm for solving quadratic minimization problems
// subject to box bounds. Solves x = argmin g'(x-xn) + 1/2 (x-xn)' H (x-xn).
void BoundQuadCD(double *g, double *H, double *x, int num_par, double mu, 
                 int maxit, double *xmin, double *xmax)
{
   int      i, j, it, found = 0;
   double   xn2, err, dj, xj;
   double   hes, tmp;
   double   etol = 1.0e-9;

   it = 0; 
   while ((found == 0) && (it < maxit)) {
      xn2 = vecnorm2(x, num_par);
      err = 0.0;
      for (j = 0; j < num_par; j++) {
         // Estimate the step size
         tmp = mu;
         hes = H[j + j * num_par] + tmp;
         if (fabs(hes) > 0.0)
            dj = g[j] / hes;
         else
            dj = 0.0;
         xj = x[j] - dj;

         // Apply bounds to the estimate
         if (xj < xmin[j])
            xj = xmin[j];
         else if (xj > xmax[j])
            xj = xmax[j];
         dj = xj - x[j];
         x[j] = xj;

         // Update gradient
         for (i = 0; i < num_par; i++)
            g[i] += H[i + j * num_par] * dj;
         g[j] += tmp * dj;

         // Update the error
         err += dj * dj;
      }
      ++it;

      // Check for convergence
      if (sqrt(err) <= etol * (xn2 + etol))
         found = 1;
   }
}

//------------------------------------------------------------------------------
// lema_gsn
//------------------------------------------------------------------------------
// The Levenberg-Marquardt algorithm to minimize the least squares problem
// subject to linear bounds. This version uses a gradient descent approach.
void lema_gsn(double *w, double *y, double *f, int num_y, double *p, int num_p, 
              void *param, void (*func)(double *, void *, double *),
              void (*jacf)(double *, void *, double *, int *, double *),
              double *plb, double *pub, int *psens, int maxit)
{
   int      i, j, k, it, n;
   double   mu, v, tau, rho, maxh;
   int      num_par;
   double   *pnew, *pest, *pdif, *plb_sens, *pub_sens;
   double   *J, *g, *H, *gnew;
   double   etol = 1e-9;
   double   F, Fnew, L;
   int      subit = 100;
   double   tmp, h1, h2;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (j = 0; j < num_p; j++)
      if (psens[j] == 1) ++num_par;

   // Memory allocation
   pdif = (double*) malloc(sizeof(double) * num_par);
   pnew = (double*) malloc(sizeof(double) * num_par);
   pest = (double*) malloc(sizeof(double) * num_par);
   plb_sens = (double*) malloc(sizeof(double) * num_par);
   pub_sens = (double*) malloc(sizeof(double) * num_par);
   J = (double*) malloc(sizeof(double) * num_y * num_par);
   g = (double*) malloc(sizeof(double) * num_par);
   H = (double*) malloc(sizeof(double) * num_par * num_par);
   gnew = (double*) malloc(sizeof(double) * num_par);

   // Initialize the parameter estimates
   getkin(p, psens, num_par, pest);
   getkin(plb, psens, num_par, plb_sens);
   getkin(pub, psens, num_par, pub_sens);
   (*jacf)(p, param, f, psens, J);

   // Compute the initial cost function value
   F = 0.0;
   for (i = 0; i < num_y; i++)
      F += w[i] * (y[i] - f[i]) * (y[i] - f[i]) * 0.5;

   // Initialize search parameters
   v = 2.0;
   tau = 1.0e-3;
   maxh = 0;
   for (j = 0; j < num_par; j++) {
      tmp = vecnorm2(J + num_y * j, num_y);
      tmp *= tmp;
      if (maxh < tmp) maxh = tmp;
   }
   mu = tau * maxh;

   // Iterative loop for optimization
   it = 0; n = 0;
   rho = 1.0;
   while (it < maxit) {
      // Calculate gradient and approximate Hessian
      if (rho > 0.0) {
         for (j = 0; j < num_par; j++) {
            g[j] = 0.0;
            for (i = 0; i < num_y; i++) {
               h1 = w[i] * (f[i] - y[i]);
               g[j] +=  h1 * J[i + j * num_y];
            }
            for (k = 0; k < num_par; k++) {
               H[j + k * num_par] = 0.0;
               for (i = 0; i < num_y; i++) {
                  h2 = w[i];
                  H[j + k * num_par] += h2 * J[i + j * num_y] * J[i + k * num_y];
               }
            }
         }
      }

      // Estimate the parameters
      for (j = 0; j < num_par; j++) {
         pnew[j] = pest[j];
         gnew[j] = g[j];
      }
      BoundQuadCD(gnew, H, pnew, num_par, mu, subit, plb_sens, pub_sens);
      for (i = 0; i < num_par; i++)
         pdif[i] = pnew[i] - pest[i];
      L = 0.0;
      for (j = 0; j < num_par; j++) {
         L += g[j] * pdif[j];
         for (i = 0; i < num_par; i++)
            L += H[j + i * num_par] * pdif[j] * pdif[i] * 0.5;
      }

      if (vecnorm2(pdif, num_par) <= etol * (vecnorm2(pest, num_par) + etol)) break;

      // Calculate the ratio rho
      setkin(pnew, num_par, psens, p);
      (*func)(p, param, f);
      Fnew = 0.0;
      for (i = 0; i < num_y; i++)
         Fnew += w[i] * (y[i] - f[i]) * (y[i] - f[i]) * 0.5;
      rho = (Fnew - F) / L;

      // Update mu and parameter estimates
      if (rho > 0) {
         for (j = 0; j < num_par; j++)
            pest[j] = pnew[j];
         (*jacf)(p, param, f, psens, J);
         F = Fnew;
         ++it;

         tmp = 2 * rho - 1;
         tmp = 1 - tmp * tmp * tmp;
         if (tmp < 1.0 / 3.0) //must be 1.0/3.0, 1/3 doesn't work
            tmp = 1.0 / 3.0;
         mu *= tmp;
         v = 2.0;
      } 
      else {
         mu *= v;
         v *= 2.0;
      }

      ++n;
      if (n > 1000) {
         printf("stopped: maximum iteration number exceeds 1000!\n");
         setkin(pest, num_par, psens, p);
         (*func)(p, param, f);
         break;
      }
   }

   // Free allocated memory
   free(pdif);
   free(pnew);
   free(pest);
   free(plb_sens);
   free(pub_sens);
   free(J);
   free(g);
   free(H);
   free(gnew);
}

//------------------------------------------------------------------------------
// kconv_liver_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve using a two-tissue kinetic model with 
// liver-specific parameters.
void kconv_liver_tac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, ka, fa;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t;
   int      num_time;

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_pv = (double*) malloc(sizeof(double) * num_time);
   cp   = (double*) malloc(sizeof(double) * num_time);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);
   c_t  = (double*) malloc(sizeof(double) * num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Transform kinetic parameters into exponential parameters
      vb = p[0+j*7];
      k1 = p[1+j*7];
      k2 = p[2+j*7];
      k3 = p[3+j*7];
      k4 = p[4+j*7];
      ka = p[5+j*7];
      fa = p[6+j*7];
      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;

      // Calculate dispersion and input concentration
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca, num_time, td, c_pv);
      for (i = 0; i < num_time; i++) {
         cp[i] = (1 - fa) * c_pv[i] + fa * ca[i];
      }

      // Calculate the exponential components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // Calculate the concentrations
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;
      for (i = 0; i < num_time; i++) {
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * cp[i];
      }

      // Average the values over frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(cp);
   free(c_pv);
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
}

//------------------------------------------------------------------------------
// kconv_liver_jac
//------------------------------------------------------------------------------
// Calculates the time activity curves and sensitivity functions for the 
// two-tissue compartment model specific to liver kinetic modeling.
void kconv_liver_jac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, ka, fa;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 7; i++) {
      if (psens[i] == 1)
         ++num_par;
   }

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_pv = (double*) malloc(sizeof(double) * num_time);
   cp   = (double*) malloc(sizeof(double) * num_time);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);
   c_t  = (double*) malloc(sizeof(double) * num_time);
   s_t  = (double*) malloc(sizeof(double) * num_time * num_par);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Transform kinetic parameters into exponential parameters
      vb = p[0 + j * 7];
      k1 = p[1 + j * 7];
      k2 = p[2 + j * 7];
      k3 = p[3 + j * 7];
      k4 = p[4 + j * 7];
      ka = p[5 + j * 7];
      fa = p[6 + j * 7];
      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;

      // Calculate dispersion
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca, num_time, td, c_pv);
      for (i = 0; i < num_time; i++) {
         cp[i] = (1 - fa) * c_pv[i] + fa * ca[i];
      }

      // Calculate the exponential components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // Calculate concentrations and sensitivity functions
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;
      for (i = 0; i < num_time; i++) {
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * cp[i];
      }

      // Average the time activity curve over the frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);

      // Calculate sensitivity functions for each parameter
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + cp[i];
         s_t += num_time;
      }
      if (psens[1] == 1 || psens[2] == 1) { 
         f1 = 1 / d * (k4 + k3 - a1);
         f2 = 1 / d * (a2 - k4 - k3);
      }
      if (psens[1] == 1) { // Sensitivity wrt k1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[2] == 1 || psens[3] == 1) { 
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a2);
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[3] == 1 || psens[4] == 1) { 
         f1 = 1 / d * (a1 + a2 - k3 - k4);
         f2 = -f1;
      }
      if (psens[3] == 1) { // Sensitivity wrt k3
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[4] == 1) { // Sensitivity wrt k4
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a2);
         for (i = 0; i < num_time; i++)
            s_t[i] = -(1 - vb) * (f1 * c_a1[i] + f2 * c_a2[i]);
         s_t += num_time;
      }
      if (psens[5] == 1 || psens[6] == 1) { 
         for (i = 0; i < num_time; i++) {
            cp[i] = ca[i] - c_pv[i];
         }
         f1 = k1 / d * (k4 - a1);
         f2 = k1 / d * (a2 - k4);
      }
      if (psens[5] == 1) { // Sensitivity wrt ka
         tmp = ka + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_t);
         for (i = 0; i < num_time; i++) {
            c_t[i] = (1.0 - fa) * c_t[i];
         }
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_t, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_t, num_time, td, c_a2);
         for (i = 0; i < num_time; i++) {
            s_t[i] = (1 - vb) * ((f1 + b1) * c_a1[i] + (f2 + b2) * c_a2[i]) + vb * c_t[i];
         }
         s_t += num_time;
      }
      if (psens[6] == 1) { // Sensitivity wrt fa
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a2);
         for (i = 0; i < num_time; i++) {
            s_t[i] = (1 - vb) * ((f1 + b1) * c_a1[i] + (f2 + b2) * c_a2[i]) + vb * cp[i];
         }
         s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Average the sensitivity functions over the frame duration
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }

   // Free allocated memory
   free(cp);
   free(c_pv);
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------




