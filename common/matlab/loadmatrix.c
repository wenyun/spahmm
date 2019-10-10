#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "io.h"

void readRow(char* line,
             double* Sub,
             int nRows,
             int nCols,
             int row,
             double* Cols) {

  int i = 0;
  int c = 0;
  char* tok;
  int flag = 1;

  while(1) {
    if (flag) {
      tok = strtok(line, " \t");
      flag = 0;
    } else {
      tok = strtok(NULL, " \t");
    }
    
    if (i+1 == (int) Cols[c]) {
      Sub[row + c * nRows] = atof(tok);
      if (++c == nCols) {
        break;
      }
    }
    ++i;
    if(tok == NULL || *tok == '\n') 
      break;
  }

  if (c != nCols) {
    mexErrMsgTxt("column index out of range\n");
  }
}

bool checkIncreasingOrder(double* Array, int N) {
  
  int i;
  for (i = 1; i < N; i++) {
    if (Array[i] < Array[i-1]) {
      return false;
    }
  }
  return true;
}

/* *  
 * This function loads sub matrix by specified indices
 * @ param (file, rows, cols)
 * file - the file containing whole matrix
 * rows - a matrix containing all rows for submatrix
 * cols - a matrix containing all cols for submatrix
 * The index must be in increasing order (TODO change it)
 * */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (nrhs < 3 || nrhs > 3) {
    mexErrMsgTxt("only supports two parameters\n");
  }
  
  max_line_len = 1024;
  line = malloc(max_line_len * sizeof(char));

  char* MatFile = mxArrayToString(prhs[0]);
  double* Rows = mxGetPr(prhs[1]);
  int nRows = mxGetM(prhs[1]) * mxGetN(prhs[1]);
  double* Cols = mxGetPr(prhs[2]);
  int nCols = mxGetM(prhs[2]) * mxGetN(prhs[2]);

  if (nRows == 0 || nCols == 0) {
    return;
  }

  plhs[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
  double* Sub = mxGetPr(plhs[0]);

  if (!checkIncreasingOrder(Rows, nRows) || 
      !checkIncreasingOrder(Cols, nCols)) {
    mexErrMsgTxt("Index is not in increasing order\n");
  }

  FILE* fp = fopen((char*) MatFile, "r");
  int i = 0;
  int r = 0;
  line = readline(fp);
  while (line != NULL) {      
    if (i+1 == (int) Rows[r]){
      readRow(line, Sub, nRows, nCols, r, Cols);
      if (++r == nRows) {
        break;
      }
    }
    ++i;
    line = readline(fp);
  }

  if (r != nRows) {
    mexErrMsgTxt("row index out of range\n");
  }

  fclose(fp);
  free(line);
  mxFree(MatFile);
}


