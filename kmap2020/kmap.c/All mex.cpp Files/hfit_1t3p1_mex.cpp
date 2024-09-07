#include "mex.h"
#include "kinlib.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Usage: kfit_2t(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td, tac2, w2)      
//          
{ 
  int i, j;
  double *tac, *w, *scant, *cp, *wb, *tac2, *w2;
  double dk;
  double *pinit;
  int num_frm, num_vox, num_par, np, nw, nw2;
  int psens[5];
  double *temp;
  int maxit; 
  double *p, *c, *cj, *wj, *pj, *cfit, *wj2, *cj2;
  double *plb, *pub;
  double td;
  KMODEL_T km;
  
  // input   
  tac = mxGetPr(prhs[0]);  
  num_frm = (int) mxGetM(prhs[0]);
  num_vox = (int) mxGetN(prhs[0]);
  w = mxGetPr(prhs[1]);
  nw= mxGetN(prhs[1]);   
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
  tac2 = mxGetPr(prhs[12]); 
  w2 = mxGetPr(prhs[13]);   
  nw2= mxGetN(prhs[13]);
  
  // model parameters
  km.dk = dk;
  km.td = td;
  km.cp = cp;
  km.wb = wb;
  km.num_frm = num_frm;
  km.num_vox = 1;
  km.scant = scant;
  km.tacfunc = kconv_1t3p1_tac;
  km.jacfunc = kconv_1t3p1_jac;

  // lable sensitive parameters
  if (num_par==1)
    printf("WARNING: pinit should be a column vector or a matrix!");
  for (i=0; i<num_par; i++) {
    psens[i] = (int) temp[i];
  }
  
  // memory for output
  plhs[0] = mxCreateDoubleMatrix(num_par, num_vox, mxREAL);
  p = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(num_frm, num_vox, mxREAL);
  c = mxGetPr(plhs[1]);
    
  // initialize parameter p
  if (np==1) {
    for (i=0; i<num_par; i++)
      for (j=0; j<num_vox; j++)
        p[i+j*num_par] = pinit[i];
  } else if (np==num_vox) {
    for (i=0; i<num_par; i++)
      for (j=0; j<num_vox; j++)
        p[i+j*num_par] = pinit[i+j*num_par];
  }

  // voxel-wise fitting
  for (j=0; j<num_vox; j++) {
    
    cj = tac + j*num_frm;
    if (nw==num_vox)
      wj = w + j*num_frm;
    else
      wj = w;    
    pj = p + j*num_par;
    cj2 = tac2 + j*num_frm;
    if (nw2==num_vox)
      wj2 = w2 + j*num_frm;
    else
      wj2 = w2; 
      
    cfit = c + j*num_frm;
    lema_png(wj, cj, cfit, wj2, cj2, num_frm, pj, num_par, &km, tac_eval, 
             jac_eval, plb, pub, psens, maxit);
  }
}
