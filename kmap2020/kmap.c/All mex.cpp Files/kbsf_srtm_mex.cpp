#include "mex.h"
#include "kinlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Usage: u = kbsf_srtm_mex(p, c, scant, dk, td)      
//         
{ 
   double  *p, *c, *u;
   double  *scant;
   int      num_par, num_vox, num_frm, num_time;
   double   dk, td;
  
   // input   
   p = mxGetPr(prhs[0]);  
   num_par = (int) mxGetM(prhs[0]);
   num_vox = (int) mxGetN(prhs[0]);
   c = mxGetPr(prhs[1]);  
   num_frm = (int) mxGetM(prhs[1]);
   scant = mxGetPr(prhs[2]);      
   dk = mxGetScalar(prhs[3]);
   td = mxGetScalar(prhs[4]);
    
   // check
   if (num_vox!=(int) mxGetN(prhs[1]))
      printf("unmatched size of kinetic parameters and tac");
   if (num_par==1)
      printf("WARNING: pinit should be a column vector or a matrix!");
  
   // memory for output
   num_time = (int) (scant[2*num_frm-1]/td);
   plhs[0] = mxCreateDoubleMatrix(num_time, 1, mxREAL);
   u = mxGetPr(plhs[0]);
  
   // the blood spread function  
   kconv_srtm_bsf(p, dk, scant, td, c, num_frm, num_vox, u);   

}
