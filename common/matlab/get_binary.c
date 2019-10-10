#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "binary.h"

int checkRange(double* Array, int M, int N, int Max) {

  int i, j;
  for (i = 0; i < M; i++) {
    for(j = 0; j < N; j++) {
      if (Array[i + j * M] > Max) {
        mexErrMsgTxt("Index out of range");
        return 0;
      }
    }
  }
  return 1;
}

unsigned int getValue(unsigned char* BinaryArray, int M, int N, int i, int j) {
  
  int Index = j / BIT_PER_BYTE;
  int Offset = j % BIT_PER_BYTE;
  unsigned int Value = ((BinaryArray[Index * M + i] & bitmask[Offset]) > 0);

  return Value;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (nrhs != 3) {
    mexErrMsgTxt("only support 3 parameters");
  }

  unsigned char* DataMat = mxGetData(prhs[0]);
  int M = mxGetM(prhs[0]);
  int N = mxGetN(prhs[0]);

  double* I = mxGetPr(prhs[1]);
  int IM = mxGetM(prhs[1]);
  int IN = mxGetN(prhs[1]);
  double* J = mxGetPr(prhs[2]);
  int JM = mxGetM(prhs[2]);
  int JN = mxGetN(prhs[2]);

  checkRange(I, IM, IN, M);
  checkRange(J, JM, JN, BIT_PER_BYTE * N);
  
  int SM = IM * IN;
  int SN = JM * JN;

  plhs[0] = mxCreateNumericMatrix(SM, SN, mxUINT8_CLASS, mxREAL);
  unsigned char* SubMat = mxGetData(plhs[0]);

  int i, j;
  int in, im, jn, jm;
  for (in = 0; in < IN; in++) {
    for (im = 0; im < IM; im++) {
      i = in * IM + im;
      for (jn = 0; jn < JN; jn++) {
        for (jm = 0; jm < JM; jm++) {
          j = jn * JM + jm;
          SubMat[j * SM + i] = getValue(DataMat, M, N, ((int)I[i]) - 1, ((int)J[j]) - 1);
        }
      }
    }
  }
}

