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
// kconv_2tcm_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the two-tissue compartment model.
// This model uses six parameters: vb, k1, k2, k3, k4, and t_delay.
void kconv_2tcm_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t;
   double  *cp_delay, *wb_delay;  
   int      num_time;

   // Memory allocation for intermediate calculations
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);
   cp_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   for (j=0; j<num_vox; j++) {
      // Parameter transformations
      vb = p[0+j*6];
      k1 = p[1+j*6];
      k2 = p[2+j*6];
      k3 = p[3+j*6];
      k4 = p[4+j*6];
      t_delay = p[5+j*6];
      k234 = k2+k3+k4;
      d = sqrt(k234*k234 - 4*k2*k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;

      // get delayed cp/wb
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(cp, num_time, t_delay, td, cp_delay);
      // Convolution with exponential functions
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a2);

      // Compute tissue concentration
      if (d==0) d = 1e9;
      f1 = k1/d * (k4 - a1);
      f2 = k1/d * (a2 - k4);
      b1 = k1/d * k3;
      b2 = -b1;

      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * wb_delay[i];
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
   free(cp_delay);
   free(wb_delay);
}

//------------------------------------------------------------------------------
// kconv_2tcm_jac
//------------------------------------------------------------------------------
// Calculate the time activity curves and the sensitivity functions 
// (Jacobian) for the two-tissue compartment model.
void kconv_2tcm_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, k3, k4, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   double  *c_a1_0, *c_a2_0, *c_f_0, *c_b_0, *c_t_0, *wb_delay;       
   double  *grad_t_delay;                          
   int      num_time;
   int      num_par;
  
   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 6; i++) {
      if (psens[i] == 1)
         ++num_par;
   }
  
   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a1 = (double*) malloc(sizeof(double) * num_time);
   c_a2 = (double*) malloc(sizeof(double) * num_time);
   c_f  = (double*) malloc(sizeof(double) * num_time);
   c_b  = (double*) malloc(sizeof(double) * num_time);
   c_a1_0 = (double*) malloc(sizeof(double)*num_time);
   c_a2_0 = (double*) malloc(sizeof(double)*num_time);
   c_f_0  = (double*) malloc(sizeof(double)*num_time);
   c_b_0  = (double*) malloc(sizeof(double)*num_time);
   c_t_0  = (double*) malloc(sizeof(double)*num_time);
   wb_delay  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double) * num_time);
   s_t  = (double*) malloc(sizeof(double) * num_time * num_par);
	grad_t_delay  = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {
	
      // Transform kinetic parameters into exponential parameters
      vb = p[0 + j * 6];
      k1 = p[1 + j * 6];
      k2 = p[2 + j * 6];
      k3 = p[3 + j * 6];
      k4 = p[4 + j * 6];
      t_delay = p[5+j*6];

      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;
	
      // exp components without time delay
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1_0);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2_0);
      // c_a1, c_a2, wb with time delay
      time_delay_tac(c_a1_0, num_time, t_delay, td, c_a1);
      time_delay_tac(c_a2_0, num_time, t_delay, td, c_a2);
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);

      // Calculate concentrations
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;	
      for (i = 0; i < num_time; i++) {
         c_f_0[i] = f1*c_a1_0[i] + f2*c_a2_0[i];
         c_b_0[i] = b1*c_a1_0[i] + b2*c_a2_0[i];
         c_t_0[i] = (1-vb) * (c_f_0[i] + c_b_0[i]) + vb*wb[i];
         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * wb_delay[i];
      }

      // Sensitivity functions (Jacobian) computation
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + wb_delay[i];
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
      if (psens[5] == 1) { // wrt time delay
         time_delay_jac(c_t_0, num_time, t_delay, td, grad_t_delay);
	      for (i=0; i<num_time; i++){ 
            s_t[i] = - grad_t_delay[i];
         } 
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
   free(c_a1_0);
   free(c_a2_0);
   free(c_f_0);
   free(c_b_0);
   free(wb_delay);
   free(c_t_0);
   free(grad_t_delay);
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
   double   vb, R1, k2, BP, t_delay;
   double   k2a, a;
   double  *c_a, *c_t;
   double  *cr0_delay, *wb_delay;
   int      num_time;

   // Memory allocation for intermediate calculations
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);
   cr0_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   for (j = 0; j < num_vox; j++) {
      // Parameter transformations
      vb = p[0 + j * 5];
      R1 = p[1 + j * 5];
      k2 = p[2 + j * 5];
      BP = p[3 + j * 5];
      t_delay = p[4+j*5];

      // get delayed cr0/wb
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(cr0, num_time, t_delay, td, cr0_delay);
      // Convolution with exponential functions
      k2a = k2 / (1.0 + BP);
      a = k2a + dk;
      kconv_exp(1.0, a, cr0_delay, num_time, td, c_a);

      // Compute tissue concentration
      for (i = 0; i < num_time; i++)
         c_t[i] = (1 - vb) * (R1 * cr0_delay[i] + (k2 - R1 * k2a) * c_a[i]) + vb * wb_delay[i];

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
   free(cr0_delay);
   free(wb_delay);
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
   double   vb, R1, k2, BP, t_delay;
   double   k2a;
   double   *c_a, *c_b, *c_t, *s_t;
   double   *c_a_0, *c_b_0, *c_t_0, *wb_delay, *cr0_delay;
   double   *grad_t_delay;
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
   c_a = (double*) malloc(sizeof(double) * num_time);
   c_b = (double*) malloc(sizeof(double) * num_time);
   c_a_0 = (double*) malloc(sizeof(double) * num_time);
   c_b_0 = (double*) malloc(sizeof(double) * num_time);
   c_t_0 = (double*) malloc(sizeof(double) * num_time);
   wb_delay = (double*) malloc(sizeof(double) * num_time);
   cr0_delay = (double*) malloc(sizeof(double) * num_time);

   c_t = (double*) malloc(sizeof(double) * num_time);
   s_t = (double*) malloc(sizeof(double) * num_time * num_par);

   grad_t_delay  = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Assign parameters
      vb = p[0 + j * 5];
      R1 = p[1 + j * 5];
      k2 = p[2 + j * 5];
      BP = p[3 + j * 5];
      t_delay = p[4 + j * 5];

      // Compute k2a
      k2a = k2 / (1.0 + BP);
      
      // Calculate the exponential components
      kconv_exp(1.0, k2a + dk, cr0, num_time, td, c_a_0);

      // shift c_a, wb, cr0 for time delay
      time_delay_tac(c_a_0, num_time, t_delay, td, c_a);
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(cr0, num_time, t_delay, td, cr0_delay);

      // Calculate concentration c_t
      for (i = 0; i < num_time; i++) {
         c_b_0[i] = R1 * cr0[i] + (k2 - R1 * k2a) * c_a_0[i];
         c_t_0[i] = (1 - vb) * c_b_0[i] + vb * wb[i];

         c_b[i] = R1 * cr0_delay[i] + (k2 - R1 * k2a) * c_a[i];
         c_t[i] = (1 - vb) * c_b[i] + vb * wb_delay[i];
      }

      // Average the values over frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);

      // Calculate sensitivity functions for each parameter
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = -c_b[i] + wb_delay[i];
         s_t += num_time;
      }
      if (psens[1] == 1) { // Sensitivity wrt R1
         for (i = 0; i < num_time; i++)
            s_t[i] = (1 - vb) * (cr0_delay[i] - k2a * c_a[i]);
         s_t += num_time;
      }
      if (psens[2] == 1) { // Sensitivity wrt k2
         for (i = 0; i < num_time; i++)
            c_t[i] = cr0_delay[i] - c_b[i] / (1.0 + BP);
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
      if (psens[3] == 1) { // wrt time delay
         time_delay_jac(c_t_0, num_time, t_delay, td, grad_t_delay);
	      for (i=0; i<num_time; i++){ 
            s_t[i] = - grad_t_delay[i];
         } 
	      s_t += num_time;
      }
      s_t -= num_time * num_par;

      // Average the sensitivity functions over frame duration
      frame(scant, td, s_t, num_frm, num_par, st + j * num_frm * num_par);
   }

   // Free allocated memory
   free(c_a);
   free(c_b);
   free(c_a_0);
   free(c_b_0);
   free(c_t_0);
   free(cr0_delay);
   free(wb_delay);
   free(grad_t_delay);

   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------
// kconv_1tcm_tac
//------------------------------------------------------------------------------
// Calculates the time activity curve (TAC) using the one-tissue compartment model.
// This model uses four parameters: vb, k1, k2, and t_delay.
void kconv_1tcm_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
{
   int      i, j;
   double   vb, k1, k2, t_delay;
   double   a;
   double  *c_a, *c_t;
   double  *cp_delay, *wb_delay;      
   int      num_time;

   // Memory allocation
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);

   cp_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   for (j = 0; j < num_vox; j++) {
      // Assign parameters
      vb = p[0 + j * 4];
      k1 = p[1 + j * 4];
      k2 = p[2 + j * 4];
      t_delay = p[3+j*4];

      time_delay_tac(cp, num_time, t_delay, td, cp_delay);
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      // Convolution with exponential functions
      a = k2 + dk;
      kconv_exp(k1, a, cp_delay, num_time, td, c_a);

      // Compute tissue concentration
      for (i = 0; i < num_time; i++)
         c_t[i] = (1-vb)*c_a[i] + vb*wb_delay[i];

      // Frame-averaged activity
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);
   }

   // Free allocated memory
   free(c_a);
   free(c_t);
   free(cp_delay);
   free(wb_delay);
}

//------------------------------------------------------------------------------
// kconv_1tcm_jac
//------------------------------------------------------------------------------
// Calculates the time activity curves and the sensitivity functions for the 
// one-tissue kinetic model.
void kconv_1tcm_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
{
   int      i, j;
   double   vb, k1, k2, t_delay;
   double   a;
   double  *c_a, *c_t, *s_t;
   double  *c_a_0, *c_t_0, *wb_delay;  
   double  *grad_t_delay;                         
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 4; i++) {
      if (psens[i] == 1) ++num_par;
   }

   // Allocate memory for intermediate calculations
   num_time = (int) (scant[2*num_frm-1] / td);
   c_a = (double*) malloc(sizeof(double) * num_time);
   c_t = (double*) malloc(sizeof(double) * num_time);
   s_t = (double*) malloc(sizeof(double) * num_time * num_par);

   c_a_0 = (double*) malloc(sizeof(double)*num_time);
   c_t_0 = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);
   grad_t_delay  = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Assign parameters
      vb = p[0 + j * 4];
      k1 = p[1 + j * 4];
      k2 = p[2 + j * 4];
      t_delay = p[3 + j * 4];

      // Calculate the exponential component
      a = k2 + dk;
      kconv_exp(k1, a, cp, num_time, td, c_a_0);
      // shift c_a, wb for time delay
      time_delay_tac(c_a_0, num_time, t_delay, td, c_a);
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      // Calculate the time activity curve
      for (i = 0; i < num_time; i++){
         c_t_0[i] = (1-vb) * c_a_0[i] + vb * wb[i];
         c_t[i] = (1 - vb) * c_a[i] + vb * wb_delay[i];
      }

      // Calculate sensitivity functions
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = -c_a[i] + wb_delay[i];
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
      if (psens[3] == 1) { // wrt time delay
         time_delay_jac(c_t_0, num_time, t_delay, td, grad_t_delay);
	      for (i=0; i<num_time; i++){ 
            s_t[i] = - grad_t_delay[i];
         } 
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
   free(c_a_0);
   free(c_t_0);
   free(wb_delay);
   free(grad_t_delay);
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
   double   vb, k1, k2, k3, k4, ka, fa, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t;
   double   *ca_delay, *wb_delay;              
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
   ca_delay = (double*) malloc(sizeof(double)*num_time);
   wb_delay = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Transform kinetic parameters into exponential parameters
      vb = p[0+j*8];
      k1 = p[1+j*8];
      k2 = p[2+j*8];
      k3 = p[3+j*8];
      k4 = p[4+j*8];
      ka = p[5+j*8];
      fa = p[6+j*8];
      t_delay = p[7 + j*8];

      k234 = k2 + k3 + k4;
      d = sqrt(k234 * k234 - 4 * k2 * k4);
      a1 = (k234 - d) / 2;
      a2 = (k234 + d) / 2;
      // get delayed cp/wb
      time_delay_tac(wb, num_time, t_delay, td, wb_delay);
      time_delay_tac(ca, num_time, t_delay, td, ca_delay);

      // Calculate dispersion and input concentration
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca_delay, num_time, td, c_pv);
      for (i = 0; i < num_time; i++) {
         cp[i] = (1 - fa) * c_pv[i] + fa * ca_delay[i];
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
   free(ca_delay);
   free(wb_delay);
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
   double   vb, k1, k2, k3, k4, ka, fa, t_delay;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   double   *cp_delay, *c_f_0, *c_b_0, *c_t_0, *c_a1_0, *c_a2_0, *grad_t_delay;
   int      num_time;
   int      num_par;

   // Determine the number of sensitive parameters
   num_par = 0;
   for (i = 0; i < 8; i++) {
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

   cp_delay   = (double*) malloc(sizeof(double) * num_time);
   c_a1_0 = (double*) malloc(sizeof(double) * num_time);
   c_a2_0 = (double*) malloc(sizeof(double) * num_time);
   c_f_0  = (double*) malloc(sizeof(double) * num_time);
   c_b_0  = (double*) malloc(sizeof(double) * num_time);
   c_t_0  = (double*) malloc(sizeof(double) * num_time);

   c_t  = (double*) malloc(sizeof(double) * num_time);
   s_t  = (double*) malloc(sizeof(double) * num_time * num_par);

   grad_t_delay  = (double*) malloc(sizeof(double)*num_time);

   // Iterate over each voxel
   for (j = 0; j < num_vox; j++) {

      // Transform kinetic parameters into exponential parameters
      vb = p[0 + j * 8];
      k1 = p[1 + j * 8];
      k2 = p[2 + j * 8];
      k3 = p[3 + j * 8];
      k4 = p[4 + j * 8];
      ka = p[5 + j * 8];
      fa = p[6 + j * 8];
      t_delay = p[7 + j * 8];

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
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1_0);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2_0);
      // shift c_a1, c_a2, cp for time delay
      time_delay_tac(c_a1_0, num_time, t_delay, td, c_a1);
      time_delay_tac(c_a2_0, num_time, t_delay, td, c_a2);
      time_delay_tac(cp, num_time, t_delay, td, cp_delay);

      // Calculate concentrations and sensitivity functions
      if (d == 0) d = 1e9;
      f1 = k1 / d * (k4 - a1);
      f2 = k1 / d * (a2 - k4);
      b1 = k1 / d * k3;
      b2 = -b1;
      for (i = 0; i < num_time; i++) {
         c_f_0[i] = f1 * c_a1_0[i] + f2 * c_a2_0[i];
         c_b_0[i] = b1 * c_a1_0[i] + b2 * c_a2_0[i];
         c_t_0[i] = (1 - vb) * (c_f_0[i] + c_b_0[i]) + vb * cp[i];

         c_f[i] = f1 * c_a1[i] + f2 * c_a2[i];
         c_b[i] = b1 * c_a1[i] + b2 * c_a2[i];
         c_t[i] = (1 - vb) * (c_f[i] + c_b[i]) + vb * cp_delay[i];
      }

      // Average the time activity curve over the frame duration
      frame(scant, td, c_t, num_frm, 1, ct + j * num_frm);

      // Calculate sensitivity functions for each parameter
      if (psens[0] == 1) { // Sensitivity wrt vb
         for (i = 0; i < num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + cp_delay[i];
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
         time_delay_tac(cp, num_time, t_delay, td, cp_delay);
         f1 = k1 / d * (k4 - a1);
         f2 = k1 / d * (a2 - k4);
      }
      if (psens[5] == 1) { // Sensitivity wrt ka
         tmp = ka + dk;
         kconv_exp(1.0, tmp, cp_delay, num_time, td, c_t);
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
         kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, cp_delay, num_time, td, c_a2);
         for (i = 0; i < num_time; i++) {
            s_t[i] = (1 - vb) * ((f1 + b1) * c_a1[i] + (f2 + b2) * c_a2[i]) + vb * cp_delay[i];
         }
         s_t += num_time;
      }
      if (psens[7] == 1) { // wrt time delay
         time_delay_jac(c_t_0, num_time, t_delay, td, grad_t_delay);
	      for (i=0; i<num_time; i++){ 
            s_t[i] = - grad_t_delay[i];
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

   free(cp_delay);
   free(c_f_0);
   free(c_b_0);
   free(c_t_0);
   free(c_a1_0);
   free(c_a2_0);
   free(grad_t_delay);

   free(c_t);
   free(s_t);
}

//------------------------------------------------------------------------------

