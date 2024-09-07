#include "mex.h"
#include "kinlib.h"


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Usage: ktac_2t(par, scant, blood, wblood, dk, td)      
//          
{   
  double *par, *scant, *b, *wb;
  double dk, td;
  int num_par, num_vox, num_frm; 
  double *c, *s;
  int psens[3] = {1, 1, 1};
  
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
  c = mxGetPr(plhs[0]);
  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(num_frm, num_par, mxREAL);
    s = mxGetPr(plhs[1]);
  }
    
  // voxel-wise 
  if (nlhs==1)
    kconv_1t3p_tac(par, dk, scant, td, b, wb, num_frm, num_vox, c);
  else if (nlhs>1)
    kconv_1t3p_jac(par, dk, scant, td, b, wb, num_frm, num_vox, c, psens, s);
}
