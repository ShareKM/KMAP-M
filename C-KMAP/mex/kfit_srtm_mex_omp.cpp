#include "mex.h"
#include "kinlib.h"
#include <omp.h>

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This file implements the fitting of a Simplified Reference Tissue Model (SRTM)
// using the Levenberg-Marquardt algorithm with OpenMP for parallel processing 
// within the MATLAB environment.
//
// Usage: 
// kfit_srtm_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
//
// Compilation Instruction:
// mex kfit_srtm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp 
// CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
//
// This will produce a MEX file named 'kfit_srtm_mex_omp', which you can call from MATLAB 
// as kfit_srtm_mex_omp(...) with the same arguments as described above.
//
// Input parameters:
// - tac: Time activity curve (TAC) data.
// - w: Weights for the TAC data.
// - scant: Scan time data.
// - blood: Blood data.
// - wblood: Whole blood data.
// - dk: Decay constant.
// - pinit: Initial parameters for the model.
// - lb: Lower bounds for the parameters.
// - ub: Upper bounds for the parameters.
// - psens: Sensitivity matrix for the parameters.
// - maxit: Maximum number of iterations for the fitting algorithm.
// - td: Time duration for the scan.
//
// Output:
// - p: Estimated parameters.
// - c: Fitted curve.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int i, j;
    double *tac, *w, *scant, *cp, *wb, *w1;
    double dk;
    double *pinit;
    int num_frm, num_vox, num_par, np, nw;
    int psens[5];
    double *temp;
    int maxit; 
    double *p, *c;
    double *plb, *pub;
    double td;
    KMODEL_T km;

    // Retrieve input arguments from MATLAB
    tac = mxGetPr(prhs[0]);  
    num_frm = (int) mxGetM(prhs[0]);
    num_vox = (int) mxGetN(prhs[0]);
    nw = mxGetN(prhs[1]);  
    w1 = mxGetPr(prhs[1]);  
    scant = mxGetPr(prhs[2]);      
    cp = mxGetPr(prhs[3]); 
    wb = mxGetPr(prhs[4]); 
    dk = mxGetScalar(prhs[5]);
    pinit = mxGetPr(prhs[6]); 
    num_par = mxGetM(prhs[6]);
    np = (int) mxGetN(prhs[6]);
    plb = mxGetPr(prhs[7]);
    pub = mxGetPr(prhs[8]);
    temp = mxGetPr(prhs[9]);
    maxit = (int) mxGetScalar(prhs[10]);
    td = mxGetScalar(prhs[11]);

    // Set up the kinetic model parameters
    km.dk = dk;
    km.td = td;
    km.cp = cp;
    km.wb = wb;
    km.num_frm = num_frm;
    km.num_vox = 1; // Process one voxel at a time
    km.scant = scant;
    km.tacfunc = kconv_srtm_tac; // TAC function for SRTM model
    km.jacfunc = kconv_srtm_jac; // Jacobian function for SRTM model

    // Label sensitive parameters
    if (num_par == 1) {
        mexWarnMsgTxt("WARNING: pinit should be a column vector or a matrix!");
    }
    for (i = 0; i < num_par; i++) {
        psens[i] = (int)temp[i];
    }

    // Allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(num_par, num_vox, mxREAL);
    p = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(num_frm, num_vox, mxREAL);
    c = mxGetPr(plhs[1]);

    // Initialize parameter p
    if (np == 1) {
        for (i = 0; i < num_par; i++) {
            for (j = 0; j < num_vox; j++) {
                p[i + j * num_par] = pinit[i];
            }
        }
    } 
    else if (np == num_vox) {
        for (i = 0; i < num_par; i++) {
            for (j = 0; j < num_vox; j++) {
                p[i + j * num_par] = pinit[i + j * num_par];
            }
        }
    }

    // Prepare weights
    if (nw == 1) {
        w = (double*) malloc(sizeof(double) * num_frm * num_vox);
        for (i = 0; i < num_frm; i++) {
            for (j = 0; j < num_vox; j++) {
                w[i + j * num_frm] = w1[i];    
            }
        }
    } 
    else {
        w = w1;
    }

    // Display the total number of threads available
    i = omp_get_max_threads();
    mexPrintf("Total number of threads: %d\n", i);

    #pragma omp parallel
    {
        int tid, m;
        double *cj, *wj, *pj, *cfit;
        
        // Allocate memory for each thread
        cj = (double*) malloc(sizeof(double) * num_frm);
        wj = (double*) malloc(sizeof(double) * num_frm);
        pj = (double*) malloc(sizeof(double) * num_par);
        cfit = (double*) malloc(sizeof(double) * num_frm);

        tid = omp_get_thread_num();
        
        #pragma omp for nowait
        // Voxel-wise fitting
        for (j = 0; j < num_vox; j++) {

            // Copy TAC and weights for the current voxel
            for (m = 0; m < num_frm; m++) {
                cj[m] = tac[m + j * num_frm];
                wj[m] = w[m + j * num_frm];
            }

            // Copy initial parameters for the current voxel
            for (m = 0; m < num_par; m++) {
                pj[m] = p[m + j * num_par];
            }

            // Perform Levenberg-Marquardt fitting
            kmap_levmar(cj, wj, num_frm, pj, num_par, &km, tac_eval, jac_eval, plb, pub, 
                        psens, maxit, cfit);

            // Copy fitted curve and parameters back to the output
            for (m = 0; m < num_frm; m++) {
                c[m + j * num_frm] = cfit[m];
            }
            for (m = 0; m < num_par; m++) {
                p[m + j * num_par] = pj[m];
            }
        }

        // Free memory allocated for each thread
        if (cj) free(cj);
        if (wj) free(wj);
        if (pj) free(pj);
        if (cfit) free(cfit);
    }

    // Free the weight array if it was dynamically allocated
    if ((w) && (nw == 1)) {
        free(w);
    }
}
