#include "mex.h"
#include "kinlib.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          
{   
  double *par, *scant, *b, *wb;
  double dk, td;
  int num_par, num_vox, num_frm; 
  double *ct, *cf, *cb;
  int psens[7] = {1, 1, 1, 1, 1, 1, 1};
  
  // input   
  par = mxGetPr(prhs[0]);  
  num_par = (int) mxGetM(prhs[0]);
  num_vox = (int) mxGetN(prhs[0]);
  scant = mxGetPr(prhs[1]);      
  num_frm = (int) mxGetM(prhs[1]);
  b = mxGetPr(prhs[2]); 
  wb = mxGetPr(prhs[3]); 
  dk = mxGetScalar(prhs[4]);
  td = mxGetScalar(prhs[5]);
  
  // output
  plhs[0] = mxCreateDoubleMatrix(num_frm, num_vox, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(num_frm, num_vox, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(num_frm, num_vox, mxREAL);
  ct = mxGetPr(plhs[0]);
  cf = mxGetPr(plhs[1]);  
  cb = mxGetPr(plhs[2]);

  // voxel-wise 
  kconv_liver_mcd(par, dk, scant, td, b, wb, num_frm, num_vox, ct, cf, cb);
  
}
