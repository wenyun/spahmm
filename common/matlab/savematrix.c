#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (nrhs < 3 || nrhs > 3) {
    mexErrMsgTxt("only supports three parameters\n");
  }

  double* Matrix = mxGetPr(prhs[0]);
  int M = mxGetM(prhs[0]);
  int N = mxGetN(prhs[0]);
  
  if (M == 0 || N == 0) {
    return;
  }

  char* MatFile = mxArrayToString(prhs[1]);
  char* Format = mxArrayToString(prhs[2]);
  char* PlaceHolder = malloc(10 * sizeof(char));
  if (strcmp(Format, " %d") == 0) {
    mxFree(Format);
    Format = malloc(10 * sizeof(char));
    strcpy(Format, " %.0f");
  } else {    
    strcpy(PlaceHolder, Format);
    mxFree(Format);
    Format = PlaceHolder;
  }

  FILE* fp = fopen((char*) MatFile, "w");
  int i, j;
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      fprintf(fp, Format, Matrix[i + j * M]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  mxFree(MatFile);
  free(Format);
}


 
