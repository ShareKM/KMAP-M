#include "mex.h"
#include <stdio.h>
#include <omp.h>

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
// kfit_test_omp.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
{
   int maxth; 
   
   maxth = omp_get_max_threads();
   printf("Number of threads: %d\n", maxth);
   
   #pragma omp parallel
   {
   int count;
   int th_id;
   th_id = omp_get_thread_num();
   printf("Hello World from master thread %d\n", th_id);
   
   #pragma omp for
   for( count = 0; count < 6; count++ )
   {
      th_id = omp_get_thread_num();
      printf("Hello World from thread %d, %d\n", th_id, count);
   }
   }
   printf( "Finished\n" );
}
