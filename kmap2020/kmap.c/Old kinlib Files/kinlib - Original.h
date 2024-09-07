#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// #include <matrix.h>

/*----------------------------------------------------------------------------*/
/*                              Kinetic Modeling                              */
//------------------------------------------------------------------------------
typedef struct {
   double  dk;      // decay constant
   double  td;      // time interval for convolution calculation (in sec)
   double *cp;      // concentration in plasma (in sec)
   double *wb;      // concentration in whole blood (in sec)
   int     num_frm; // number of time frame;
   int     num_vox; // number of voxels;
   double *scant;   // scan time = [t_start, t_end] (in sec)
   void  (*tacfunc)(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct); 
                    // time acitivity curve
   void  (*jacfunc)(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);	
                    // Jacobian (sensitivity functions)
} KMODEL_T;
//------------------------------------------------------------------------------
void tac_eval(double *p, void *param, double *tac);
//------------------------------------------------------------------------------
void jac_eval(double *p, void *param, double *tac, int *psens, double *jac);
//------------------------------------------------------------------------------
void frame(double *scant, double td, double *c_t, int num_frm, int num_c, double *c);
//------------------------------------------------------------------------------
void kconv_exp(double k1, double k2, double *u, int num_time, double td, double *c);
//------------------------------------------------------------------------------
void kconv_2t5p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_2t5p_MCD(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, double *cf, double *cb);
//------------------------------------------------------------------------------
void kconv_2t5p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);      				
//------------------------------------------------------------------------------
void kconv_2t5pv_tac(double *p, double dk, double *scant, double td, double *cp, 
                     double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_2t5pv_jac(double *p, double dk, double *scant, double td, double *cp, 
                     double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                     double *st); 
//------------------------------------------------------------------------------
void kconv_2t5pk_tac(double *p, double dk, double *scant, double td, double *cp, 
                     double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_2t5pk_jac(double *p, double dk, double *scant, double td, double *cp, 
                     double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                     double *st); 
//------------------------------------------------------------------------------
void kconv_2t6p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_2t6p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st); 
//------------------------------------------------------------------------------
void kconv_srtm_tac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_srtm_jac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);
//------------------------------------------------------------------------------
void kconv_1t3p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_1t3p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void kconv_1t3p1_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_1t3p1_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);
//------------------------------------------------------------------------------
void kconv_1t4p_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_1t4p_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);
//------------------------------------------------------------------------------

/*----------------------------------------------------------------------------*/
/*                         Optimization   Algorithm                           */
//------------------------------------------------------------------------------
void kmap_levmar(double *y, double *w, int num_frm, double *pinit, int num_p,
                 void *param, void (*func)(double *, void *, double *),
                 void (*jacf)(double *, void *, double *, int *, double *),
                 double *plb, double *pub, int *psens, int maxit, double *ct);
//------------------------------------------------------------------------------
void boundpls_cd(double *y, double *w, double *a, double alpha, int num_y,
                 int num_x, double *xlb, double *xub, int maxit, double *x, 
                 double *r);
//------------------------------------------------------------------------------
void setkin(double *x, int num_x, int *xsens, double *x0);
//------------------------------------------------------------------------------
void getkin(double *x0, int *xsens, int num_x, double *x);
//------------------------------------------------------------------------------
double vecnorm2(double *x, int num);
//------------------------------------------------------------------------------
double vecnormw(double *w, double *x, int num);
//------------------------------------------------------------------------------
void lema_gsn(double *w, double *y, double *f, int num_y, double *p, int num_p, 
              void *param, void (*func)(double *, void *, double *),
              void (*jacf)(double *, void *, double *, int *, double *),
              double *plb, double *pub, int *psens, int maxit);
//------------------------------------------------------------------------------
void lema_psn(double *w, double *y, double *f, int num_y, double *p, int num_p, 
              void *param, void (*func)(double *, void *, double *),
              void (*jacf)(double *, void *, double *, int *, double *),
              double *plb, double *pub, int *psens, int maxit);
//------------------------------------------------------------------------------
void lema_png(double *w, double *y, double *f, double *wg, double *yg, int num_y, 
              double *p, int num_p, 
              void *param, void (*func)(double *, void *, double *),
              void (*jacf)(double *, void *, double *, int *, double *),
              double *plb, double *pub, int *psens, int maxit);
//------------------------------------------------------------------------------             
void BoundQuadCD(double *g, double *H, double *x, int num_par, double mu, 
                 int maxit, double *xmin, double *xmax);
//------------------------------------------------------------------------------       
//------------------------------------------------------------------------------
void kconv_liver_tac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------  
void kconv_liver_mcd(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, double *cf, double *cb);
//------------------------------------------------------------------------------      
void kconv_liver_jac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);
//------------------------------------------------------------------------------    
void kconv_dbif_tac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------   
void kconv_dbif_jac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);
//------------------------------------------------------------------------------
void kconv_2tss_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);
//------------------------------------------------------------------------------
void kconv_2tss_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);    
//------------------------------------------------------------------------------
