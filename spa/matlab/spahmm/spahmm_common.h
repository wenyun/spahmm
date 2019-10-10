#pragma once

#include "numeric.h"
#include "binary.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

unsigned int getAllele(const unsigned char* BinaryArray, int M, int N, int i, int j) {
  
  int Index = j / BIT_PER_BYTE;
  int Offset = j % BIT_PER_BYTE;
  unsigned int Value = ((BinaryArray[Index * M + i] & bitmask[Offset]) > 0);

  return Value;
}

void getParameter(const mxArray* param,
                  int* nHap,
                  int* nSnp,
                  double* lambda,
                  double* omega,
                  double* rho, 
                  double** recombProb,
                  double** nonRecombProb,
                  int* nSample) {
  
  int i;
  const char* fieldName;
  const mxArray *fieldPtr;
  double* ptr;
  int nFields;
  
  nFields = mxGetNumberOfFields(param);
  for (i = 0; i < nFields; i++) {
    fieldName = mxGetFieldNameByNumber(param, i);
    fieldPtr = mxGetFieldByNumber(param, 0, i);
    ptr = mxGetPr(fieldPtr);

    if (strcmp(fieldName, "nSnp") == 0) {
      *nSnp = (int) ptr[0];
    } else if (strcmp(fieldName, "nHap") == 0) {
      *nHap = (int) ptr[0];
    } else if (strcmp(fieldName, "lambda") == 0) {
      *lambda = ptr[0];
    } else if (strcmp(fieldName, "omega") == 0) {
      *omega = ptr[0];
    } else if (strcmp(fieldName, "rho") == 0) {
      *rho = ptr[0];
    } else if (strcmp(fieldName, "recombProb") == 0) {
      *recombProb = ptr;
    } else if (strcmp(fieldName, "nonRecombProb") == 0) {
      *nonRecombProb = ptr;
    } else if (strcmp(fieldName, "nSample") == 0) {
      *nSample = (int) ptr[0];
    }   
  }
}

void getHapref(const mxArray* hapref,
               unsigned char** haps,
               double** loc,
               int* nHap,
               int* nSnp) {
  
  int i;
  const char* fieldName;
  const mxArray *fieldPtr;
  unsigned char* ptr;
  int nFields;
  
  nFields = mxGetNumberOfFields(hapref);
  for (i = 0; i < nFields; i++) {
    fieldName = mxGetFieldNameByNumber(hapref, i);
    fieldPtr = mxGetFieldByNumber(hapref, 0, i);

    if (strcmp(fieldName, "haps") == 0) {
      *haps = mxGetData(fieldPtr);
    } else if (strcmp(fieldName, "loc") == 0) {
      *loc = mxGetPr(fieldPtr);
    } else if (strcmp(fieldName, "nHap") == 0) {
      *nHap = (int) mxGetScalar(fieldPtr);
    } else if (strcmp(fieldName, "nSnp") == 0) {
      *nSnp = (int) mxGetScalar(fieldPtr);
    }
  }
}

void compute_distance(double* dist, 
                      const double* loc, 
                      const double* refloc, 
                      const double* index, 
                      const int nRefHap,
                      const int nHap) {

  int i;
  for (i = 0; i < nHap; i++) {    
    dist[i] = sqrt(pow(loc[0] - refloc[(int)index[i] - 1 + 0], 2.0) + 
                   pow(loc[1] - refloc[(int)index[i] - 1 + nRefHap], 2.0));
  }
}

void compute_transition(double* tranProbBase, 
                        const double lambda,
                        const double* dist,
                        const int nHap) {
  
  int i;
  double logZ = logsum_array_weight(dist, nHap, -lambda);
  for (i = 0; i < nHap; i++) {
    tranProbBase[i] = - lambda * dist[i] - logZ;
  }
}
