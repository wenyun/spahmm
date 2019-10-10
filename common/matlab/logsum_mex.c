#include <float.h>
#include <math.h>
#include "mex.h"
#include "numeric.h"

/*
 * logsum_mex(A, B) 
 * 
 * This function performs the logsum operation for various input
 * logsum_mex(A), will do logsum for entire array
 * logsum_mex(A, b), will do logsum(a, b) for every element a
 * logsum_mex(A, B), will do elementwise logsum
 * */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int i;
  if (nrhs == 1) {
    int M = mxGetM(prhs[0]);
    int N = mxGetN(prhs[0]);
    if (M > 1 && N > 1) {
      mexErrMsgTxt("X is not an array\n");
    }

    double* logX = mxGetPr(prhs[0]);
    double LogSum = -DBL_MAX;    
    for (i = 0; i < M*N; i++) {
      LogSum = logsum_xy(LogSum, logX[i]);
    }

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* Result = mxGetPr(plhs[0]);
    Result[0] = LogSum;

  } else if (nrhs == 2) {
    int M = mxGetM(prhs[0]);
    int N = mxGetN(prhs[0]);
    int M2 = mxGetM(prhs[1]);
    int N2 = mxGetN(prhs[1]);

    double* logX = mxGetPr(prhs[0]);
    
    if (M2 == 1 && N2 == 1) {

      double logA = mxGetScalar(prhs[1]);

      plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
      double* Result = mxGetPr(plhs[0]);
      for (i = 0; i < M*N; i++) {
        Result[i] = logsum_xy(logX[i], logA);
      }
    } else if (M2 == M && N2 == N) {
      double* logY = mxGetPr(prhs[1]);
      
      plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
      double* Result = mxGetPr(plhs[0]);
      for (i = 0; i < M*N; i++) {
        Result[i] = logsum_xy(logX[i], logY[i]);
      }
    } else {
      mexErrMsgTxt("The inputs are not feasible\n");
    }
  }
}

