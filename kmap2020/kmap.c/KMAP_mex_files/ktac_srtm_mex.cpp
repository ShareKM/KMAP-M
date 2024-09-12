#include "mex.h"
#include "kinlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This file implements the computation of the time activity curve (TAC) and its
// Jacobian for a simplified reference tissue model (SRTM) within the MATLAB 
// environment.
//
// Usage:
// ktac_srtm(par, scant, blood, wblood, dk, td)
//
// Compilation Instruction:
// mex ktac_srtm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_srtm
//
// This will produce a MEX file named 'ktac_srtm', which you can call from MATLAB 
// as ktac_srtm(...) with the same arguments as described above.
//
// Input parameters:
// - par: Model parameters.
// - scant: Scan time data.
// - blood: Blood data.
// - wblood: Whole blood data.
// - dk: Decay constant.
// - td: Time duration for the scan.
//
// Output:
// - c: Computed time activity curve (TAC).
// - s: Jacobian matrix (if requested).
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {   

    double *par, *scant, *b, *wb;
    double dk, td;
    int num_par, num_vox, num_frm; 
    double *c, *s;
    int psens[5] = {1, 1, 1, 1, 1};  // Default sensitivity for 4 parameters
  
    // Retrieve input arguments from MATLAB
    par = mxGetPr(prhs[0]);  
    num_par = (int)mxGetM(prhs[0]);
    num_vox = (int)mxGetN(prhs[0]);
    scant = mxGetPr(prhs[1]);      
    num_frm = (int)mxGetM(prhs[1]);
    b = mxGetPr(prhs[2]); 
    wb = mxGetPr(prhs[3]); 
    dk = mxGetScalar(prhs[4]);
    td = mxGetScalar(prhs[5]);
  
    // Allocate memory for output
    plhs[0] = mxCreateDoubleMatrix(num_frm, num_vox, mxREAL);
    c = mxGetPr(plhs[0]);
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(num_frm, num_par, mxREAL);
        s = mxGetPr(plhs[1]);
    }
    
    // Perform voxel-wise computation
    if (nlhs == 1) {
        kconv_srtm_tac(par, dk, scant, td, b, wb, num_frm, num_vox, c);
    } else if (nlhs > 1) {
        kconv_srtm_jac(par, dk, scant, td, b, wb, num_frm, num_vox, c, psens, s);
    }
}
