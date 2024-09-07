/*
 * kinlib_models.cpp
 * 
 * This file contains all the kinetic modeling functions used to calculate time activity curves (TAC) and sensitivity functions for different kinetic models.
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

