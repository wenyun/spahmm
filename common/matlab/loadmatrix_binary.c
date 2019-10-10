#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "binary.h"
#include "io.h"

void count_file(FILE *fp, int *rows, int *columns) {
  int rw = 0;
  int cl = 0;
  int k;
  char *tok;

  while(readline(fp) != NULL) {
    k = 0;
    if (line[0] != '#') {
      tok = strtok(line, " \t");
      ++k;
      while(1) {
        tok = strtok(NULL, " \t");
        /* check '\n' as ' ' may be after the last element */
        if(tok == NULL || *tok == '\n') 
          break;
        ++k;
      }
      
      if(k == 0){
        mexErrMsgTxt("empty row\n");
      }
      if(cl == 0) {
        cl = k;
      } else {
        if(cl != k) {
          mexErrMsgTxt("different number of columns\n");
        }
      }
      rw++;
    }
  }
  
  (*rows) = rw;
  (*columns) = cl;
}

void readMatrix(FILE*fp,
                int nRows,
                int nCols,
                unsigned char* DataMat,
                int M,
                int N) {

  int i, j;
  unsigned char buf;
  int tok;

  for (i = 0; i < nRows; i++) {
    readline(fp);
    
    buf = 0;
    for (j = 0; j < nCols; j++) {
      if (j == 0) {
        tok = atoi(strtok(line, " \t"));
      } else {
        tok = atoi(strtok(NULL, " \t"));
      }

      if (tok == 1) {
        buf |= bitmask[j % BIT_PER_BYTE];
      }

      if ((j + 1) % BIT_PER_BYTE == 0) {
        DataMat[((j + 1) / BIT_PER_BYTE - 1) * M + i] = buf;
        buf = 0;
      }
    }

    /* some residual byte */
    if (nCols < BIT_PER_BYTE * N) {
      DataMat[(N-1) * M + i] = buf;
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int M;
  int N;
  int nRows;
  int nCols; 

  if (nrhs != 1) {
    mexErrMsgTxt("only support one parameters\n");
  }
  if (nlhs != 3) {
    mexErrMsgTxt("must be three output variables\n");
  }
  
  max_line_len = 1024;
  line = malloc(max_line_len * sizeof(char));

  char* MatFile = mxArrayToString(prhs[0]);
  FILE* fp = fopen(MatFile, "r");
  count_file(fp, &nRows, &nCols);
  rewind(fp);

  M = nRows;
  N = (nCols - 1) / BIT_PER_BYTE + 1;
   
  plhs[0] = mxCreateNumericMatrix(M, N, mxUINT8_CLASS, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
   
  unsigned char* DataMat = mxGetData(plhs[0]);
  double* Rows = mxGetPr(plhs[1]);
  double* Cols = mxGetPr(plhs[2]);
  
  *Rows = (double) nRows;
  *Cols = (double) nCols;
 
  readMatrix(fp, nRows, nCols, DataMat, M, N); 

  fclose(fp);
  free(line);
  mxFree(MatFile);
}

