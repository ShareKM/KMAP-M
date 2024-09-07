#include "mex.h"
#include "kinlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This file implements the fitting of a two-tissue kinetic model (2T5P) using
// the Levenberg-Marquardt algorithm within the MATLAB environment.
//
// Usage: 
// kfit_2t5p(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td)
//
// Compilation Instruction:
// mex kfit_2t5p_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2t5p
//
// This will produce a MEX file named 'kfit_2t5p', which you can call from MATLAB 
// as kfit_2t5p(...) with the same arguments as described above.
//
// Input parameters:
// - tac: Time activity curve (TAC) data.
// - w: Weights for the TAC data.
// - scant: Scan time data.
// - blood: Blood data.
// - wblood: Whole blood data.
// - dk: Decay constant.
// - pinit: Initial parameters for the model.
// - plb: Lower bounds for the parameters.
// - pub: Upper bounds for the parameters.
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
    double *tac, *w, *scant, *cp, *wb;
    double dk;
    double *pinit;
    int num_frm, num_vox, num_par, np, nw;
    int psens[5];
    double *temp;
    int maxit; 
    double *p, *c, *cj, *wj, *pj, *cfit;
    double *plb, *pub;
    double td;
    KMODEL_T km;

    // Retrieve input arguments from MATLAB
    tac = mxGetPr(prhs[0]);  
    num_frm = (int) mxGetM(prhs[0]);
    num_vox = (int) mxGetN(prhs[0]);
    w = mxGetPr(prhs[1]);
    nw = mxGetN(prhs[1]);   
    scant = mxGetPr(prhs[2]);      
    cp = mxGetPr(prhs[3]); 
    wb = mxGetPr(prhs[4]); 
    dk = mxGetScalar(prhs[5]);
    pinit = mxGetPr(prhs[6]); 
    num_par = (int) mxGetM(prhs[6]);
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
    km.tacfunc = kconv_2t5p_tac; // TAC function for 2T5P model
    km.jacfunc = kconv_2t5p_jac; // Jacobian function for 2T5P model

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

    // Voxel-wise fitting
    for (j = 0; j < num_vox; j++) {
    
        cj = tac + j * num_frm;
        if (nw == num_vox) {
            wj = w + j * num_frm;
        } else {
            wj = w;    
        }
        pj = p + j * num_par;
    
        cfit = c + j * num_frm;
        lema_gsn(wj, cj, cfit, num_frm, pj, num_par, &km, tac_eval, jac_eval, plb, pub, psens, maxit);
    }
}
