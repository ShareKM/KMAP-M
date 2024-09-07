#include "kinlib.h"
#include "mex.h"

//------------------------------------------------------------------------------
void tac_eval(double *p, void *param, double *tac)
//------------------------------------------------------------------------------
{ 
   KMODEL_T *par;
   par = (KMODEL_T *) param;
   (*(par->tacfunc))(p, par->dk, par->scant, par->td, par->cp, par->wb, 
                     par->num_frm, par->num_vox, tac);
}

//------------------------------------------------------------------------------
void jac_eval(double *p, void *param, double *tac, int *psens, double *jac)
//------------------------------------------------------------------------------
{
   KMODEL_T *par;
   par = (KMODEL_T *) param;
   (*(par->jacfunc))(p, par->dk, par->scant, par->td, par->cp, par->wb, 
                     par->num_frm, par->num_vox, tac, psens, jac);
}

//------------------------------------------------------------------------------
void kconv_2t5p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
//------------------------------------------------------------------------------
// calculate the time acitivity curve using the two-tissue kinetic model. 
{
   int      i, j;
   double   vb, k1, k2, k3, k4;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t;
   int      num_time;
	
   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);
	
   for (j=0; j<num_vox; j++) {
	
      // transform k's to exp parameters
      vb = p[0+j*5];
      k1 = p[1+j*5];
      k2 = p[2+j*5];
      k3 = p[3+j*5];
      k4 = p[4+j*5];
      k234 = k2+k3+k4;
      d = sqrt(k234*k234-4*k2*k4);
      a1 = (k234-d)/2;
      a2 = (k234+d)/2;

      // exp components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);
	
      // [c_f, c_b] and c
      if (d==0) d = 1e9;
      f1 = k1/d*(k4-a1);
      f2 = k1/d*(a2-k4);
      b1 = k1/d*k3;
      b2 = -b1;
      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1-vb) * (c_f[i] + c_b[i]) + vb*wb[i];
      }
	
      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
   }
  	
   // release memory
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
	
}


//------------------------------------------------------------------------------
void kconv_2t5p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
//------------------------------------------------------------------------------
// calculate the time activity curves in second and the sensitivity functions 
//	for the two-tissue compartment model.
{
   int      i, j;
   double   vb, k1, k2, k3, k4;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double  *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   int      num_time;
   int      num_par;
  
   // determine the number of sensitive parameters
   num_par = 0;
   for (i=0; i<5; i++) {
      if (psens[i]==1)
         ++num_par;
   }
  
   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);
   s_t  = (double*) malloc(sizeof(double)*num_time*num_par);
	
   // voxel-wise
   for (j=0; j<num_vox; j++) {
	
      // transform k's to exp parameters
      vb = p[0+j*5];
      k1 = p[1+j*5];
      k2 = p[2+j*5];
      k3 = p[3+j*5];
      k4 = p[4+j*5];
      k234 = k2+k3+k4;
      d = sqrt(k234*k234-4*k2*k4);
      a1 = (k234-d)/2;
      a2 = (k234+d)/2;
	
		// exp components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);
	
      // [c_f, c_b] and c
      if (d==0) d = 1e9;
      f1 = k1/d*(k4-a1);
      f2 = k1/d*(a2-k4);
      b1 = k1/d*k3;
      b2 = -b1;	
      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1-vb) * (c_f[i] + c_b[i]) + vb*wb[i];
      }
		 	
      // calculate the sensitivity functions
      if (psens[0]==1) { // wrt vb
         for (i=0; i<num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + wb[i];
         s_t += num_time;
      }
      if (psens[1]==1 || psens[2]==1) { 
         f1 = 1/d*(k4+k3-a1);
         f2 = 1/d*(a2-k4-k3);
      }
      if (psens[1]==1) { // wrt k1
         for (i=0; i<num_time; i++)
            s_t[i] = (1-vb) * ( f1*c_a1[i] + f2*c_a2[i]);
         s_t += num_time;
      }
      if (psens[2]==1 || psens[3]==1) { 
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a2);
      }
      if (psens[2]==1) { // wrt k2
         for (i=0; i<num_time; i++)
            s_t[i] = -(1-vb) * ( f1*c_a1[i] + f2*c_a2[i] );
         s_t += num_time;
      }
      if (psens[3]==1 || psens[4]==1) { 
         f1 = 1/d*(a1+a2-k3-k4);
         f2 = -f1;
      }
		if (psens[3]==1) { // wrt k3
			for (i=0; i<num_time; i++)
				s_t[i] = (1-vb) * ( f1*c_a1[i] + f2*c_a2[i] );
			s_t += num_time;
		}
      if (psens[4]==1) { // wrt k4
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a2);
         for (i=0; i<num_time; i++)
            s_t[i] = -(1-vb) * ( f1*c_a1[i] + f2*c_a2[i] );
         s_t += num_time;
      }
      s_t -= num_time*num_par;

      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
      frame(scant, td, s_t, num_frm, num_par, st+j*num_frm*num_par);
   }
  	
   // release memory
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
   free(s_t);
}



//------------------------------------------------------------------------------
void kconv_srtm_tac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct)
//------------------------------------------------------------------------------
// calculate the time acitivity curve using the simplified reference tissue 
// model
{
   int      i, j;
   double   vb, R1, k2, BP;
   double   k2a, a;
   double  *c_a, *c_t;
   int      num_time;
        
   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);
 
   // voxel-wise
   for (j=0; j<num_vox; j++) {

      // transform k's to exp parameters
      vb = p[0+j*4];
      R1 = p[1+j*4];
      k2 = p[2+j*4];
      BP = p[3+j*4];

      // exp components
      k2a = k2/(1.0+BP);
      a = k2a + dk;
      kconv_exp(1.0, a, cr0, num_time, td, c_a);

      // c_t
      for (i=0; i<num_time; i++)
         c_t[i] = (1-vb) * ( R1 * cr0[i] + (k2-R1*k2a)*c_a[i] ) + vb*wb[i];
	   
      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
   }
   
   // release memory
   free(c_a);
   free(c_t);
	
}


//------------------------------------------------------------------------------
void kconv_srtm_jac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
//------------------------------------------------------------------------------
// calculate the time activity curves in second and the sensitivity functions
//	for the simplified reference tissue model.
{
   int      i, j;
   double   vb, R1, k2, BP;
   double   k2a;
   double  *c_a, *c_b, *c_t, *s_t;
   int      num_time;
   int      num_par;
	
   // determine the number of sensitive parameters
   num_par = 0;
   for (i=0; i<5; i++) {
      if (psens[i]==1) ++num_par;
   }

   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_b = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);
   s_t = (double*) malloc(sizeof(double)*num_time*num_par);

  	for (j=0; j<num_vox; j++) {
  	
      // transform k's to exp parameters
      vb = p[0+j*4];
      R1 = p[1+j*4];
      k2 = p[2+j*4];
      BP = p[3+j*4];
	
      // exp components
      k2a = k2/(1+BP);
      kconv_exp(1.0, k2a+dk, cr0, num_time, td, c_a);
	
      // c_t
      for (i=0; i<num_time; i++) {
         c_b[i] = R1*cr0[i] + (k2-R1*k2a)*c_a[i];
         c_t[i] = (1-vb)*c_b[i] + vb*wb[i];
      }
	   
	   // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
      
      // calculate the sensitivity functions
      if (psens[0]==1) { // wrt vb
         for (i=0; i<num_time; i++)
            s_t[i] = - c_b[i] + wb[i];
         s_t += num_time;
      }
      if (psens[1]==1) { // wrt R1
         for (i=0; i<num_time; i++)
            s_t[i] = (1-vb) * ( cr0[i] - k2a*c_a[i] );
         s_t += num_time;
      }
      if ( psens[2]==1 ) { // wrt k2
         for (i=0; i<num_time; i++)
            c_t[i] = cr0[i] - c_b[i] / (1.0+BP);
         kconv_exp(1, k2a+dk, c_t, num_time, td, c_a);
         for (i=0; i<num_time; i++)
            s_t[i] = (1-vb) * c_a[i];
         s_t += num_time;
      }
      if (psens[3]==1) { // wrt BP
         kconv_exp(k2a/(1+BP), k2a+dk, c_b, num_time, td, c_a);
         for (i=0; i<num_time; i++)
            s_t[i] = (1-vb) * c_a[i];
         s_t += num_time;
      }
      s_t -= num_time*num_par;
	
      // average in frame duration
      frame(scant, td, s_t, num_frm, num_par, st+j*num_frm*num_par);
   }
   
   // release memory
   free(c_a);
   free(c_b);
   free(c_t);
   free(s_t);
	
}

//------------------------------------------------------------------------------
void kconv_1t3p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct)
//------------------------------------------------------------------------------
// calculate the time acitivity curves using the one tissue  kinetic model.
{
   int      i, j;
   double   vb, k1, k2;
   double   a;
   double  *c_a, *c_t;
   int      num_time;
	
   // collect memory
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);
  	
   // voxel-wise
   for (j=0; j<num_vox; j++) {
  	
      // assign parameters
      vb = p[0+j*3];
      k1 = p[1+j*3];
      k2 = p[2+j*3];
	
      // exp components
      a = k2 + dk;
      kconv_exp(k1, a, cp, num_time, td, c_a);
	
      // tac
      for (i=0; i<num_time; i++)
         c_t[i] = (1-vb)*c_a[i] + vb*wb[i];
	
      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
   }
  	
   // release memory
   free(c_a);
   free(c_t);
	
}

//------------------------------------------------------------------------------
void kconv_1t3p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
//------------------------------------------------------------------------------
// calculate the time acitivity curves and the sensitivity functions for the one
// tissue kinetic model.
{
   int      i, j;
   double   vb, k1, k2;
   double   a;
   double  *c_a, *c_t, *s_t;
   int      num_time;
   int      num_par;
	
   // determine the number of sensitive parameters
   num_par = 0;
   for (i=0; i<3; i++) {
      if (psens[i]==1) ++num_par;
   }
  
   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_a = (double*) malloc(sizeof(double)*num_time);
   c_t = (double*) malloc(sizeof(double)*num_time);
   s_t = (double*) malloc(sizeof(double)*num_time*num_par);
  
   // voxel-wise
   for (j=0; j<num_vox; j++) {
  	
      // assign parameters
      vb = p[0+j*3];
      k1 = p[1+j*3];
      k2 = p[2+j*3];
	
      // exp components
      a = k2 + dk;
      kconv_exp(k1, a, cp, num_time, td, c_a);
	
      // tac
      for (i=0; i<num_time; i++)
			c_t[i] = (1-vb)*c_a[i] + vb*wb[i];
	
      // sensitivity functions
      if (psens[0]==1) { // wrt fv
         for (i=0; i<num_time; i++)
            s_t[i] = - c_a[i] + wb[i];
         s_t += num_time;
      }
      if (psens[1]==1) { // wrt k1
         for (i=0; i<num_time; i++)
            s_t[i] = (1-vb) * c_a[i]/k1;
         s_t += num_time;
      }
      if (psens[2]==1) { // wrt k2
         kconv_exp(1, a, c_a, num_time, td, s_t);
         for (i=0; i<num_time; i++)
            s_t[i] *= -(1-vb);
         s_t += num_time;
      }
      s_t -= num_time*num_par;
	
      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
      frame(scant, td, s_t, num_frm, num_par, st+j*num_frm*num_par);
   }
  	
   // release memory
   free(c_a);
   free(c_t);
   free(s_t);
	
}



//------------------------------------------------------------------------------
void frame(double *scant, double td, double *c_t, int num_frm, int num_c, double *ct)
//------------------------------------------------------------------------------
// get the average acitivity in the frame duration
{
   int      m, k;
   double   tmp, dt;
   int      num_time;
   int      j, jstart, jend;
	
   num_time = (int) (scant[2*num_frm-1]/td); 
	
   for (m=0; m<num_frm; m++) {
      jstart = (int) (scant[m]/td);
      jend = (int) (scant[m+num_frm]/td);
      dt = (double) (jend-jstart);
      for (k=0; k<num_c; k++) {
         tmp = 0.0;
         for (j = jstart; j < jend; j++) {
            tmp += c_t[j+k*num_time];
         }
         ct[m+k*num_frm] = tmp / dt;
      }
    
   }
}


//------------------------------------------------------------------------------
void kconv_exp(double k1, double k2, double *u, int num_time, double td, double *c)
//------------------------------------------------------------------------------
// compute the exponential function convoluted with the input 
//	function: c(t) = k1*exp(-k2*t) \otimes u(t). Note that t is in second.   			   
{
   int      i;
	double   ek2, prev, tmp;
	
   k1 *= td/60.0;  // e.g. td = 1 sec;
   k2 *= td/60.0;
   ek2 = exp(-k2);
   prev = 0;
   if (ek2==1.0) {
      for (i=0; i<num_time; i++) {
         prev += u[i];
         c[i] = k1 * prev;
      }
   } else {
      tmp = (1.0-ek2)/k2;
      for (i=0; i<num_time; i++) {
         prev = prev * ek2 + k1 * tmp * u[i];
         c[i] = prev;
      }
   }
}

//------------------------------------------------------------------------------
void kmap_levmar(double *y, double *w, int num_frm, double *pinit, int num_p, 
                 void *param, void (*func)(double *, void *, double *),
                 void (*jacf)(double *, void *, double *, int *, double *),
                 double *plb, double *pub, int *psens, int maxit, double *ct)
//------------------------------------------------------------------------------
// The Levenberg-Marquardt algorithm to solve the noninear least squares problem
// subject to linear bounds.
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
	
   // determine the number of sensitive parameters
   num_par = 0;
   for (j=0; j<num_p; j++)
      if (psens[j]==1) ++num_par;
  
   // collect memory space
   f = (double*) malloc(sizeof(double)*num_frm);
   fnew = (double*) malloc(sizeof(double)*num_frm);
   r = (double*) malloc(sizeof(double)*num_frm);
   h = (double*) malloc(sizeof(double)*num_par);
   pnew = (double*) malloc(sizeof(double)*num_par);
   pest = (double*) malloc(sizeof(double)*num_par);
   plb_sens = (double*) malloc(sizeof(double)*num_par);
   pub_sens = (double*) malloc(sizeof(double)*num_par);
   st = (double*) malloc(sizeof(double)*num_frm*num_par);

   // initialize the parameter estimates and their bounds 
   getkin(pinit, psens, num_par, pest);
   getkin(plb, psens, num_par, plb_sens);
   getkin(pub, psens, num_par, pub_sens);

   // initial TAC and s (the sensitivity matrix)
   (*jacf)(pinit, param, ct, psens, st);

   // the initial residual f and the cost value F
   for (i=0; i<num_frm; i++) {
      f[i] = y[i] - ct[i];
   }
   F = vecnormw(w,f,num_frm)*0.5;
  
   // initialize the parameters v, tau, mu
   v = 2.0;
   tau = 1.0e-3;
   maxaii = 0;
   for (j=0; j<num_par; j++) {
      tmp = vecnorm2(st+num_frm*j, num_frm);
      tmp *= tmp;
      if (maxaii<tmp)
         maxaii = tmp;
   }
   mu = tau * maxaii;

   // iterative loop
   it = 0; n = 0;
   while (it<maxit) {
    
      // estimate pnew and the step h
      for (j=0; j<num_par; j++)
         pnew[j] = pest[j];
      
      // coordinate least squares
      boundpls_cd(f,w,st,mu,num_frm,num_par,plb_sens,pub_sens,subit,pnew,r);
    
      // convergent or determine the new mu;
      for (i=0; i<num_par;i++) {
         h[i] = pnew[i] - pest[i];
      }
      if ( vecnorm2(h,num_par) <= etol*(vecnorm2(pest,num_par)+etol) ) 
         break; // NOTE[cannot be RETURN, otherwise free memory may not work]
		
      // calculate the ratio rho
      setkin(pnew, num_par, psens, pinit);
      (*func)(pinit, param, ct);
      for (i=0; i<num_frm; i++)
         fnew[i] = y[i] - ct[i];
      Fnew = vecnormw(w,fnew,num_frm)*0.5;
      Lnew = vecnormw(w,r,num_frm)*0.5;
      rho  = (F-Fnew) / (F-Lnew);
		//printf("it=%d, n=%d, rho=%1.2f, mu=%1.2f, F=%5.2f, Fnew=%5.2f, Lnew=%5.2f \n", it, n, rho, mu, F, Fnew, Lnew);	
      //for (j=0; j<num_par; j++)
      //   printf("%2.4f  ",pnew[j]);
      //printf("\n");
      // accept the step or not
      if ( rho > 0 ) { 
			   
         // if accept, update the estimates
         for (j=0; j<num_par; j++) {
            pest[j] = pnew[j];
         }
         (*jacf)(pinit, param, ct, psens, st);
         for (i=0; i<num_frm; i++) {
            f[i] = fnew[i];
         }
         F = Fnew;
         ++it;	
				
         // update mu
         tmp = 2*rho-1;
         tmp = 1-tmp*tmp*tmp;
         if (tmp<1.0/3.0)   // NOTE[must be 1.0/3.0, 1/3 doesn't work.]
            tmp = 1.0/3.0;				
         mu *= tmp;
         v = 2.0;
      } else {  // otherwise, find a new step
         mu *= v;
         v *= 2.0;
      }
	   ++n;
      if (n>1000) {
         printf("stopped: maximum iteration number exceeds 1000!\n");
         setkin(pest, num_par, psens, pinit);
         (*func)(pinit, param, ct);
         break;
      }
   }
	
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
void boundpls_cd(double *y, double *w, double *a, double alpha, int num_y,
                 int num_x, double *xlb, double *xub, int maxit, double *x, 
                 double *r)
//------------------------------------------------------------------------------
//	The coordinate descent algorithm to solves the penlaized least squares 
//	problem subject to box bounds: || y - A * (x-x_n) ||^2 
{
   int      i, j, ij; 
   double   err0, err, dj, xj;
   double   tmp1, tmp2;
   double   etol = 1.0e-9;
   int      it;
	
   // the residual
   for (i=0; i<num_y; i++){
      r[i] = y[i];
   }
   err0 = vecnorm2(y, num_y);
	
   // iterate
   for (it=0; it<maxit; it++) {
		
      // convergence
      err = vecnorm2(r, num_y);
      if (err/err0<etol)
         break;
		//printf("it=%d ",it);
      // over coordinate
      for (j=0; j<num_x; j++) {
			
         // 1-D optimization
         tmp1 = tmp2 = 0;
         for (i=0; i<num_y; i++) {
            ij = i+j*num_y;
            tmp1 += w[i]*a[ij]*r[i];
            tmp2 += w[i]*a[ij]*a[ij];
         }	
         tmp2 += alpha*tmp2;
         if (tmp2==0)
            dj = 0;
         else
            dj = tmp1 / tmp2;
			xj = x[j] + dj;
			//printf("%2.4f  ",x[j]);
			
         // put bounds
         if (xj<xlb[j])
            xj = xlb[j];
         else if (xj>xub[j])
            xj = xub[j];
				
         // update x and r
         dj = xj - x[j];
         x[j] = xj;
         for (i=0; i<num_y; i++)
            r[i] -= a[i+j*num_y]*dj;
      }
      //printf("\n");
   }	
}

//------------------------------------------------------------------------------
void setkin(double *x, int num_x, int *xsens, double *x0)
//------------------------------------------------------------------------------
{
   int j = 0;
   int	k = 0;
   while ( j < num_x ) {
      if ( xsens[k]==0 )
         ++k;
      else {
         x0[k] = x[j];
         ++j;
         ++k;
      }
   }
}

//------------------------------------------------------------------------------
void getkin(double *x0, int *xsens, int num_x, double *x)
//------------------------------------------------------------------------------
{
   int j = 0;
   int	k = 0;
   while ( j < num_x ) {
      if ( xsens[k]==0 )
         ++k;
      else {
         x[j] = x0[k];
         ++j;
         ++k;
      }
   }
}

//------------------------------------------------------------------------------
double vecnorm2(double *x, int num)
//------------------------------------------------------------------------------
// The norm of vector x
{
   int      i;
   double   temp = 0.0;  
  
   for (i=0; i<num; i++)
      temp += x[i]*x[i];
   return sqrt(temp);
}


//------------------------------------------------------------------------------
double vecnormw(double *w, double *x, int num)
//------------------------------------------------------------------------------
// The weighted norm of vector x
{
   int      i;
   double   temp = 0.0;  
  
   for (i=0; i<num; i++)
      temp += w[i]*x[i]*x[i];
   return temp;
}


//------------------------------------------------------------------------------
void BoundQuadCD(double *g, double *H, double *x, int num_par, double mu, 
                 int maxit, double *xmin, double *xmax)
//------------------------------------------------------------------------------
//	The coordinate descent algorithm for solving the quadratic minimization 
//	problem subject to box bounds: x = argmin g'(x-xn)+1/2(x-xn)'H(x-xn); 
{
   int      i, j, it, found;; 
   double   xn2, err, dj, xj;
   double   hes, tmp;
   double   etol = 1.0e-9;   
   
	found = 0;
   it = 0; 
   while ( (found==0) & (it<maxit) ) {
		
      xn2 = vecnorm2(x, num_par);
      err = 0.0;
      for (j=0; j<num_par; j++) {
			
			// estimate
			tmp = mu; //mu*H[j+j*num_par];  // NOTE: tmp=mu works better for direct reconstruction
			hes = H[j+j*num_par] + tmp;
			if (fabs(hes)>0.0)
			   dj = g[j] / hes;
			else
			   dj = 0.0;
			xj = x[j] - dj;
			
         // bounded estimates
         if (xj<xmin[j])
            xj = xmin[j];
         else if (xj>xmax[j])
            xj = xmax[j];
			dj = xj - x[j];
         x[j] = xj;
         
         // gradient update
         for (i=0; i<num_par; i++)
            g[i] += H[i+j*num_par] * dj;         
         g[j] += tmp * dj; 
            
         // update the difference
         err += dj*dj;
         
      }
      ++it;
      // check convergence
      if ( sqrt(err) <= etol*(xn2+etol) )
         found = 1;    
   }

}



//------------------------------------------------------------------------------
void lema_gsn(double *w, double *y, double *f, int num_y, double *p, int num_p, 
              void *param, void (*func)(double *, void *, double *),
              void (*jacf)(double *, void *, double *, int *, double *),
              double *plb, double *pub, int *psens, int maxit)
//------------------------------------------------------------------------------
// The Levenberg-Marquardt algorithm to minimize the least squares problem
// subject to linear bounds.
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
	
   // determine the number of sensitive parameters
   num_par = 0;
   for (j=0; j<num_p; j++)
      if (psens[j]==1) ++num_par;
  
   // collect memory space
   pdif = (double*) malloc(sizeof(double)*num_par);
   pnew = (double*) malloc(sizeof(double)*num_par);
   pest = (double*) malloc(sizeof(double)*num_par);
   plb_sens = (double*) malloc(sizeof(double)*num_par);
   pub_sens = (double*) malloc(sizeof(double)*num_par);
   J = (double*) malloc(sizeof(double)*num_y*num_par);
   g = (double*) malloc(sizeof(double)*num_par);
   H = (double*) malloc(sizeof(double)*num_par*num_par);
   gnew = (double*) malloc(sizeof(double)*num_par);
   
   // initialize the parameter estimates
   getkin(p, psens, num_par, pest);
   getkin(plb, psens, num_par, plb_sens);
   getkin(pub, psens, num_par, pub_sens);
   (*jacf)(p, param, f, psens, J);

   // the cost function value F
   F = 0.0;
   for (i=0; i<num_y; i++)
      F += w[i] * ( y[i] - f[i] ) * ( y[i] - f[i] ) * 0.5;
  
   // initialize the searching parameters v, tau, mu
   v = 2.0;
   tau = 1.0e-3;
   maxh = 0;
   for (j=0; j<num_par; j++) {
      tmp = vecnorm2(J+num_y*j, num_y);
      tmp *= tmp;
      if (maxh<tmp) maxh = tmp;
   }
   mu = tau * maxh;

   // iterative loop
   it = 0; n = 0;
   rho = 1.0;
   while (it<maxit) {
      
      // calculate gradient and approximate Hessian
      if (rho>0.0)
         for (j=0; j<num_par; j++) {
            g[j] = 0.0;
            for (i=0; i<num_y; i++) {
               h1 = w[i] * (f[i]-y[i]);
               g[j] +=  h1 * J[i+j*num_y];
            }
            for (k=0; k<num_par; k++) {
               H[j+k*num_par] = 0.0;
               for (i=0; i<num_y; i++) {
                  h2 = w[i];
                  H[j+k*num_par] += h2 * J[i+j*num_y] * J[i+k*num_y];
               }
            }
         }
      
      // estimate the parameters
      for (j=0; j<num_par; j++) {
         pnew[j] = pest[j];
         gnew[j] = g[j];
      }    
      BoundQuadCD(gnew, H, pnew, num_par, mu, subit, plb_sens, pub_sens);
      for (i=0; i<num_par;i++)
         pdif[i] = pnew[i] - pest[i];
      L = 0.0;
      for (j=0; j<num_par; j++) {
         L += g[j] * pdif[j];
         for (i=0; i<num_par; i++)
            L += H[j+i*num_par]*pdif[j]*pdif[i]*0.5;
      }
  
		if (vecnorm2(pdif,num_par)<=etol*(vecnorm2(pest,num_par)+etol)) break;     
                  
      // calculate the ratio rho
      setkin(pnew, num_par, psens, p);
      (*func)(p, param, f);
      Fnew = 0.0;
      for (i=0; i<num_y; i++)
         Fnew += w[i] * ( y[i] - f[i] ) * ( y[i] - f[i] ) * 0.5;;
      rho = (Fnew-F) / L;
      //printf("it=%d, n=%d, rho=%1.2f, mu=%5.2f, F=%5.2f, Fnew=%5.2f, Lnew=%5.2f \n", it, n, rho, mu, F, Fnew, L+F);
      //for (j=0; j<num_par; j++)
      //   printf("%2.4f  ",pest[j]);
      //printf("\n");
      //for (j=0; j<num_par; j++)
      //   printf("%2.4f  ",pnew[j]);
      //printf("\n");
      
      // update mu (and parameter estimates)
      if ( rho > 0 ) {  // if accept, update the estimates      
                     
         for (j=0; j<num_par; j++)
            pest[j] = pnew[j];
         (*jacf)(p, param, f, psens, J);
         F = Fnew;
         ++it;         
         
         tmp = 2*rho-1;
         tmp = 1-tmp*tmp*tmp;
         if (tmp<1.0/3.0)   // NOTE[must be 1.0/3.0, 1/3 doesn't work.]
            tmp = 1.0/3.0;				
         mu *= tmp;
         v = 2.0;
         
      } else {  // otherwise, find a new step
      
         mu *= v;
         v *= 2.0;
      }
	  
	   ++n;
      if (n>1000) {
         printf("stopped: maximum iteration number exceeds 1000!\n");
         setkin(pest, num_par, psens, p);
         (*func)(p, param, f);
         break;
      }
      
   }
	
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
void kconv_liver_tac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct)
//------------------------------------------------------------------------------
// calculate the time acitivity curve using the two-tissue kinetic model. 
{
   int      i, j;
   double   vb, k1, k2, k3, k4, ka, fa;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2,*c_f, *c_b, *c_t;
   int      num_time;
	
   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_pv = (double*) malloc(sizeof(double)*num_time);
   cp   = (double*) malloc(sizeof(double)*num_time);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);
	
   for (j=0; j<num_vox; j++) {
	
      // transform k's to exp parameters
      vb = p[0+j*7];
      k1 = p[1+j*7];
      k2 = p[2+j*7];
      k3 = p[3+j*7];
      k4 = p[4+j*7];
      ka = p[5+j*7];
      fa = p[6+j*7];
      k234 = k2+k3+k4;
      d  = sqrt(k234*k234-4*k2*k4);
      a1 = (k234-d)/2;
      a2 = (k234+d)/2;

      // dispersion
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca, num_time, td, c_pv);
      for (i=0; i<num_time; i++) {
         cp[i] = (1-fa) * c_pv[i] + fa*ca[i];
      }
      
      // exp components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

      // [c_f, c_b] and c
      if (d==0) d = 1e9;
      f1 = k1/d*(k4-a1);
      f2 = k1/d*(a2-k4);
      b1 = k1/d*k3;
      b2 = -b1;
      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1-vb) * (c_f[i] + c_b[i]) + vb*cp[i];
      }

	      
      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
   }
  	
   // release memory
   free(cp);
   free(c_pv);
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
	
}


//------------------------------------------------------------------------------
void kconv_liver_jac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st)
//------------------------------------------------------------------------------
// calculate the time activity curves in second and the sensitivity functions 
//	for the two-tissue compartment model.
{
   int      i, j;
   double   vb, k1, k2, k3, k4, ka, fa;
   double   d, a1, a2, f1, f2, b1, b2, k234, tmp;
   double   *c_pv, *cp, *c_a1, *c_a2, *c_f, *c_b, *c_t, *s_t;
   int      num_time;
   int      num_par;
  
   // determine the number of sensitive parameters
   num_par = 0;
   for (i=0; i<7; i++) {
      if (psens[i]==1)
         ++num_par;
   }
  
   // collect memory space
   num_time = (int) (scant[2*num_frm-1]/td);
   c_pv = (double*) malloc(sizeof(double)*num_time);
   cp   = (double*) malloc(sizeof(double)*num_time);
   c_a1 = (double*) malloc(sizeof(double)*num_time);
   c_a2 = (double*) malloc(sizeof(double)*num_time);
   c_f  = (double*) malloc(sizeof(double)*num_time);
   c_b  = (double*) malloc(sizeof(double)*num_time);
   c_t  = (double*) malloc(sizeof(double)*num_time);
   s_t  = (double*) malloc(sizeof(double)*num_time*num_par);
	
   // voxel-wise
   for (j=0; j<num_vox; j++) {
	
      // transform k's to exp parameters
      vb = p[0+j*7];
      k1 = p[1+j*7];
      k2 = p[2+j*7];
      k3 = p[3+j*7];
      k4 = p[4+j*7];
      ka = p[5+j*7];
      fa = p[6+j*7];
      k234 = k2+k3+k4;
      d  = sqrt(k234*k234-4*k2*k4);
      a1 = (k234-d)/2;
      a2 = (k234+d)/2;
		//printf("test ----------- %f %f %f %f %f\n", k2,k3,k4,k234,d);	
	   
	   // dispersion
      tmp = ka + dk;
      kconv_exp(ka, tmp, ca, num_time, td, c_pv);
      for (i=0; i<num_time; i++) {
         cp[i] = (1-fa) * c_pv[i] + fa*ca[i];
      }
      
		// exp components
      tmp = a1 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
      tmp = a2 + dk;
      kconv_exp(1.0, tmp, cp, num_time, td, c_a2);

	
      // [c_f, c_b] and c
      if (d==0) d = 1e9;
      f1 = k1/d*(k4-a1);
      f2 = k1/d*(a2-k4);
      b1 = k1/d*k3;
      b2 = -b1;	
      for (i=0; i<num_time; i++) {
         c_f[i] = f1*c_a1[i] + f2*c_a2[i];
         c_b[i] = b1*c_a1[i] + b2*c_a2[i];
         c_t[i] = (1-vb) * (c_f[i] + c_b[i]) + vb*cp[i];
      }
       		
      // average in frame duration
      frame(scant, td, c_t, num_frm, 1, ct+j*num_frm);
       			
      // calculate the sensitivity functions
      if (psens[0]==1) { // wrt vb
         for (i=0; i<num_time; i++)
            s_t[i] = - (c_f[i] + c_b[i]) + cp[i];
         s_t += num_time;
      }
      if (psens[1]==1 || psens[2]==1) { 
         f1 = 1/d*(k4+k3-a1);
         f2 = 1/d*(a2-k4-k3);
      }
      if (psens[1]==1) { // wrt k1
         for (i=0; i<num_time; i++)
            s_t[i] = (1-vb) * ( f1*c_a1[i] + f2*c_a2[i]);
         s_t += num_time;
      }
      if (psens[2]==1 || psens[3]==1) { 
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_f, num_time, td, c_a2);
      }
      if (psens[2]==1) { // wrt k2
         for (i=0; i<num_time; i++)
            s_t[i] = -(1-vb) * ( f1*c_a1[i] + f2*c_a2[i] );
         s_t += num_time;
      }
      if (psens[3]==1 || psens[4]==1) { 
         f1 = 1/d*(a1+a2-k3-k4);
         f2 = -f1;
      }
		if (psens[3]==1) { // wrt k3
			for (i=0; i<num_time; i++)
				s_t[i] = (1-vb) * ( f1*c_a1[i] + f2*c_a2[i] );
			s_t += num_time;
		}
      if (psens[4]==1) { // wrt k4
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_b, num_time, td, c_a2);
         for (i=0; i<num_time; i++)
            s_t[i] = -(1-vb) * ( f1*c_a1[i] + f2*c_a2[i] );
         s_t += num_time;
      }
      if (psens[5]==1 || psens[6]==1) { 
         for (i=0; i<num_time; i++) {
            cp[i] = ca[i] - c_pv[i];
         }
         f1 = k1/d*(k4-a1);
         f2 = k1/d*(a2-k4);
      }
		if (psens[5]==1) { // wrt ka
         tmp = ka + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_t);
         for (i=0; i<num_time; i++) {
            c_t[i] = (1.0-fa)*c_t[i];
         }
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, c_t, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, c_t, num_time, td, c_a2);
         for (i=0; i<num_time; i++) {
            s_t[i] = (1-vb) * ( (f1+b1)*c_a1[i] + (f2+b2)*c_a2[i] ) + vb*c_t[i];
         }
         s_t += num_time;
		}
		if (psens[6]==1) { // wrt fa
         tmp = a1 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a1);
         tmp = a2 + dk;
         kconv_exp(1.0, tmp, cp, num_time, td, c_a2);
         for (i=0; i<num_time; i++) {
            s_t[i] = (1-vb) * ( (f1+b1)*c_a1[i] + (f2+b2)*c_a2[i] ) + vb*cp[i];
         }
         s_t += num_time;
		}   
      s_t -= num_time*num_par;

      // average in frame duration
      frame(scant, td, s_t, num_frm, num_par, st+j*num_frm*num_par);
   }
  	
   // release memory
   free(cp);
   free(c_pv);
   free(c_a1);
   free(c_a2);
   free(c_f);
   free(c_b);
   free(c_t);
   free(s_t);
}

