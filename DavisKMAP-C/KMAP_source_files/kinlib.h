/*
 * kinlib.h
 * 
 * This file contains the function declarations and data structures used for kinetic modeling 
 * and optimization algorithms. It defines various functions related to the evaluation of 
 * time-activity curves (TACs), the calculation of Jacobians, and the implementation of 
 * optimization algorithms such as the Levenberg-Marquardt method.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <matrix.h>

/*----------------------------------------------------------------------------*/
/*                              Kinetic Modeling                              */
/*----------------------------------------------------------------------------*/

/*
 * KMODEL_T
 * 
 * This structure holds the parameters and function pointers required for 
 * kinetic modeling, including decay constants, concentration data, scan times, 
 * and functions for TAC and Jacobian calculations.
 */
typedef struct {
   double  dk;      // decay constant
   double  td;      // time interval for convolution calculation (in sec)
   double *cp;      // concentration in plasma (in sec)
   double *wb;      // concentration in whole blood (in sec)
   int     num_frm; // number of time frames
   int     num_vox; // number of voxels
   double *scant;   // scan time = [t_start, t_end] (in sec)
   void  (*tacfunc)(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct); 
                    // function pointer for time-activity curve calculation
   void  (*jacfunc)(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);	
                    // function pointer for Jacobian calculation (sensitivity functions)
} KMODEL_T;

/*----------------------------------------------------------------------------*/
/* Function Declarations for Kinetic Modeling */
/*----------------------------------------------------------------------------*/

/*
 * tac_eval
 * 
 * Evaluates the time-activity curve (TAC) using the given model parameters.
 */
void tac_eval(double *p, void *param, double *tac);

/*
 * jac_eval
 * 
 * Evaluates the Jacobian matrix using the given model parameters.
 */
void jac_eval(double *p, void *param, double *tac, int *psens, double *jac);

/*
 * frame
 * 
 * Computes the frame-based TAC from the scan times and input curve.
 */
void frame(double *scant, double td, double *c_t, int num_frm, int num_c, double *c);

/*
 * kconv_exp
 * 
 * Convolves an input curve with an exponential decay function.
 */
void kconv_exp(double k1, double k2, double *u, int num_time, double td, double *c);

/*
 * kconv_2tcm_tac
 * 
 * Calculates the TAC for a 2-tissue compartment model with 5 parameters.
 */
void kconv_2tcm_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);

/*
 * kconv_2tcm_jac
 * 
 * Calculates the Jacobian for a 2-tissue compartment model with 5 parameters.
 */
void kconv_2tcm_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);

/*
 * kconv_srtm_tac
 * 
 * Calculates the TAC for the Simplified Reference Tissue Model (SRTM).
 */
void kconv_srtm_tac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct);

/*
 * kconv_srtm_jac
 * 
 * Calculates the Jacobian for the Simplified Reference Tissue Model (SRTM).
 */
void kconv_srtm_jac(double *p, double dk, double *scant, double td, double *cr0, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);

/*
 * kconv_1tcm_tac
 * 
 * Calculates the TAC for a 1-tissue compartment model with 3 parameters.
 */
void kconv_1tcm_tac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct);

/*
 * kconv_1tcm_jac
 * 
 * Calculates the Jacobian for a 1-tissue compartment model with 3 parameters.
 */
void kconv_1tcm_jac(double *p, double dk, double *scant, double td, double *cp, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);

/*
 * kconv_liver_tac
 * 
 * Calculates the TAC for a liver model.
 */
void kconv_liver_tac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct);

/*
 * kconv_liver_jac
 * 
 * Calculates the Jacobian for a liver model.
 */
void kconv_liver_jac(double *p, double dk, double *scant, double td, double *ca, 
                    double *wb, int num_frm, int num_vox, double *ct, int *psens, 
                    double *st);

/*----------------------------------------------------------------------------*/
/*                         Optimization Algorithm                             */
/*----------------------------------------------------------------------------*/

/*
 * kmap_levmar
 * 
 * Implements the Levenberg-Marquardt algorithm for minimizing the least squares 
 * problem, subject to parameter bounds.
 */
void kmap_levmar(double *y, double *w, int num_frm, double *pinit, int num_p,
                 void *param, void (*func)(double *, void *, double *),
                 void (*jacf)(double *, void *, double *, int *, double *),
                 double *plb, double *pub, int *psens, int maxit, double *ct);

/*
 * boundpls_cd
 * 
 * Implements a coordinate descent algorithm for solving a quadratic minimization 
 * problem subject to box constraints.
 */
void boundpls_cd(double *y, double *w, double *a, double alpha, int num_y,
                 int num_x, double *xlb, double *xub, int maxit, double *x, 
                 double *r);

/*
 * setkin
 * 
 * Sets the values of the sensitive parameters.
 */
void setkin(double *x, int num_x, int *xsens, double *x0);

/*
 * getkin
 * 
 * Retrieves the values of the sensitive parameters.
 */
void getkin(double *x0, int *xsens, int num_x, double *x);

/*
 * vecnorm2
 * 
 * Calculates the Euclidean norm (L2 norm) of a vector.
 */
double vecnorm2(double *x, int num);

/*
 * vecnormw
 * 
 * Calculates the weighted norm of a vector.
 */
double vecnormw(double *w, double *x, int num);

/*
 * BoundQuadCD
 * 
 * Implements a coordinate descent algorithm for solving the quadratic 
 * minimization problem subject to box bounds.
 */
void BoundQuadCD(double *g, double *H, double *x, int num_par, double mu, 
                 int maxit, double *xmin, double *xmax);

/*----------------------------------------------------------------------------*/
/*                          Time   Delay  Correction                          */
/*----------------------------------------------------------------------------*/

/*
 * time_delay_jac
 * compute gradient for time delay correction
 */
void time_delay_jac(double *tac, int tac_size, double delay_time, double td, double *out);

/*
 * time_delay_tac
 * compute time-delayed TAC curve
 */
void time_delay_tac(double* input, int size, double delay_time, double td, double *out);