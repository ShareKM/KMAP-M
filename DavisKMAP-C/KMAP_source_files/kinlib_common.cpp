/*
 * kinlib_common.cpp
 * 
 * This file contains the common utility functions that are shared across different parts of the kinetic modeling and optimization process. 
 * These include functions for vector norms, setting and getting kinetic parameters, and convolution operations.
 */

#include "kinlib.h"
#include "mex.h"
#include <cmath>
#include <cstdlib>

#define min(a,b) (a < b ? a : b)
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
// time_delay_tac
//------------------------------------------------------------------------------
// compute time-delayed TAC curve
void time_delay_tac(double* tac, int tac_size, double delay_time, double td, double *out) {
/* 
   tac: time activity curve
   tac_size: number of frames
   delay_time: time delay (in seconds)
   td: time step size
   out: output tac

   Created by Yiran Wang, modified by Yansong Zhu @ UC Davis
*/
   double shiftAmount = delay_time/td;
   double fraction = shiftAmount - floor(shiftAmount);
   int shift = floor(shiftAmount);
   // Handling right shift
   if (shift >= 0) {
      int x1;
      for(int i = 0; i < min(shift, tac_size); i++) { //padding right side with 0
         x1 = i;
         out[x1] = 0;
      }
      if (shift < tac_size) { //compute values of the shifted left boundary point with linear interpolation
         out[shift] = (1-fraction) * tac[0];   
      }

      int z1,z2,z3; 
      for(int i = shift; i < tac_size-1; i++){ //compute values of the tac curve for other sample points with linear interpolation
         z1=i+1; z2=i-shift; z3=i-shift+1;
         out[z1] = (fraction) * tac[z2] + (1-fraction) * tac[z3]; 
      }
   } 
   // Handling left shift
   else {
      shift = -shift; // make shift positive for easier handling
      shift = shift - 1;

      int a1,a2,a3;
      for(int i = 0; i < tac_size - shift - 1; i++) { //compute values of the tac curve for other sample points with linear interpolation
         a1=i; a2=i+shift; a3=i+shift+1;
         out[a1] = (fraction) * tac[a2] + (1-fraction) * tac[a3];
      }
      if (shift < tac_size) { //compute values of the shifted right boundary point with linear interpolation
         int b1, b2, b3; b1=tac_size-1-shift; b2=tac_size-1; b3=tac_size-1;
         out[b1] = (fraction) * tac[b2] + (1-fraction) * tac[b3]; 
      }
        
      int c1,c2;
      for(int i = tac_size - shift; i < tac_size; i++) {
         c1=i; c2=tac_size-1;
         out[c1] = tac[c2]; // padding new points with the tail value of the original data
      }
   }
}

//------------------------------------------------------------------------------
// time_delay_jac
//------------------------------------------------------------------------------
// compute gradient for time delay correction
void time_delay_jac(double *tac, int tac_size, double delay_time, double td, double *out){
/*
   tac: time activity curve
   tac_size: number of frames
   delay_time: time delay (in seconds)
   td: time step size
   out: output gradient vector

   Created by Yiran Wang, modified by Yansong Zhu @ UC Davis
*/
int n = floor(delay_time/td);
for(int m=0; m<tac_size; m++){
   int i = m - (n + 1);
   if ((i >= 0) && (i <= (tac_size - 2))){
      out[m] = (tac[i + 1] - tac[i]) / td;
   } 
   else{
      out[m] = 0.0;    
   }
}
}
