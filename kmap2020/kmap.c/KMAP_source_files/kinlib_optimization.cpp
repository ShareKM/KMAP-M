/*
 * kinlib_optimization.cpp
 * 
 * This file contains the optimization functions used for parameter estimation in kinetic modeling. 
 * The Levenberg-Marquardt algorithm and associated helper functions are implemented here.
 */

#include "kinlib.h"
#include "mex.h"
#include <cmath>
#include <cstdlib>

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
         if (tmp < 1.0 / 3.0)
            tmp = 1.0 / 3.0;
         mu *= tmp;
         v = 2.0;
      } else {  
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
         if (tmp < 1.0 / 3.0)
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
