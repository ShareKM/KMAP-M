#include "mex.h"
#include "kinlib.h"
#include <omp.h>

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Usage: kfit_1t4p(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)  
// Compile method: mex kfit_1t4p_mex_omp.cpp kinlib.cpp kfit_test_omp.cpp 
// CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
{ 
  int		i, j;
  double	*tac, *w, *scant, *cp, *wb, *w1;
  double	dk;
  double	*pinit;
  int		num_frm, num_vox, num_par, np, nw;
  int		psens[4];
  double	*temp;
  int		maxit; 
  double	*p, *c;
  double	*plb, *pub;
  double	td;
  KMODEL_T	km;
   
  // input   
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
  
  // model parameters
  km.dk = dk;
  km.td = td;
  km.cp = cp;
  km.wb = wb;
  km.num_frm = num_frm;
  km.num_vox = 1;
  km.scant = scant;
  km.tacfunc = kconv_1t4p_tac;
  km.jacfunc = kconv_1t4p_jac;
  
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
  }
  else if (np==num_vox) {
    for (i=0; i<num_par; i++)
      for (j=0; j<num_vox; j++)
        p[i+j*num_par] = pinit[i+j*num_par];
  }

  // w
  if (nw==1) {
    w = (double*) malloc(sizeof(double)*num_frm*num_vox);
    for (i=0; i<num_frm; i++)
      for (j=0; j<num_vox; j++)
        w[i+j*num_frm] = w1[i];    
  }
  else 
    w = w1;
  
  i = omp_get_max_threads();
  printf("Total number of threads: %d\n", i);
   
  #pragma omp parallel
  {
  int    tid, m;
  double *cj, *wj, *pj, *cfit;
  cj = (double*) malloc(sizeof(double)*num_frm);
  wj = (double*) malloc(sizeof(double)*num_frm);
  pj = (double*) malloc(sizeof(double)*num_par);
  cfit = (double*) malloc(sizeof(double)*num_frm);
  
  tid = omp_get_thread_num();  
  //printf("This is thread #%d \n", tid);
    
  #pragma omp for nowait
  // voxel-wise fitting
  for (j=0; j<num_vox; j++) {
    
    for (m=0; m<num_frm; m++){
      cj[m] = tac[m+j*num_frm];
      wj[m] = w[m+j*num_frm];
    }
    
    for (m=0; m<num_par; m++)
      pj[m] = p[m+j*num_par];
    
    kmap_levmar(cj, wj, num_frm, pj, num_par, &km, tac_eval, jac_eval, plb, pub, 
                psens, maxit, cfit);
                
    for (m=0; m<num_frm; m++)
      c[m+j*num_frm] = cfit[m];
    for (m=0; m<num_par; m++)
      p[m+j*num_par] = pj[m];
  }

  if (cj) free(cj);
  if (wj) free(wj);
  if (pj) free(pj);
  if (cfit) free(cfit);  
  }
  if ((w)&&(nw==1)) free(w);

}
