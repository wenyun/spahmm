#include "mex.h"
#include "spahmm_common.h"
#include "util.h"

/* *
 * [vsum, dist, tranBase] = 
 *     spahmm_fb_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param)
 *
 * [vsum, posterior, dist, tranBase] = 
 *     spahmm_fb_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param)
 *
 * [vsum, forwardProb, backwardProb, dist, tranBase] = 
 *     spahmm_fb_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param)
 *
 * This function performs forward-backward algorithm for different task
 * */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  if (nrhs != 6) {
    mexErrMsgTxt("Incorrect number of input parameters.");
  }
  
  int i, j, k;
  
  double* loc = mxGetPr(prhs[0]);
  double* readProb = mxGetData(prhs[1]);
  unsigned char* refhap;
  int nRefHap;
  int nRefSnp;
  double* refloc;
  getHapref(prhs[2], &refhap, &refloc, &nRefHap, &nRefSnp);  
  int snpIdx = (int) mxGetScalar(prhs[3]);
  double* hapIdx = mxGetPr(prhs[4]);
  int nHap;
  int nSnp;
  int nDim = 2; /* for simplicity, using dim = 2 */
  double lambda;  
  double omega;
  double rho;
  double* recombProb;
  double* nonRecombProb;
  int nSample;
  getParameter(prhs[5],
               &nHap,
               &nSnp,
               &lambda,
               &omega,
               &rho,
               &recombProb,
               &nonRecombProb,
               &nSample);

  if (nlhs == 3) {
     /* only forward algorithm for the total likelihood */

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nHap, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nHap, 1, mxREAL);
    
    double* totalVSum = mxGetPr(plhs[0]);
    double* dist = mxGetPr(plhs[1]);
    double* tranProbBase = mxGetPr(plhs[2]);

    double* forward = malloc_double_array(nHap * 2, 0.0);
    double* emitProb = malloc_double_array(nHap, 0.0);
    double* vSum = malloc_double_array(nSnp, 0.0);
    
    compute_distance(dist, loc, refloc, hapIdx, nRefHap, nHap);
    compute_transition(tranProbBase, lambda, dist, nHap);

    for (j = 0; j < nHap; j++) {
      int refAllele = getAllele(refhap, 
                                nRefSnp,
                                nRefHap,
                                snpIdx - 1,
                                (int)hapIdx[j] - 1);
      if (refAllele == 0) {
        emitProb[j] = logsum_xy(log(1 - omega) + readProb[0],
                                log(omega) + readProb[0 + nSnp]);
      } else if (refAllele == 1) {
        emitProb[j] = logsum_xy(log(omega) + readProb[0],
                                log(1 - omega) + readProb[0 + nSnp]);
      }
      forward[j] = tranProbBase[j] + emitProb[j];
    }

    vSum[0] = logsum_array(forward, nHap);
    for (i = 1; i < nSnp; i++) {
      array_plus_constant(&forward[((i-1) % 2) * nHap], 
                          nHap,
                          nonRecombProb[snpIdx - 1 + i - 1],
                          &forward[(i % 2) * nHap]);
    
      for (j = 0; j < nHap; j++) {
        int refAllele = getAllele(refhap, 
                                  nRefSnp,
                                  nRefHap,
                                  snpIdx - 1 + i,
                                  (int)hapIdx[j] - 1); 
        if (refAllele == 0) {
          emitProb[j] = logsum_xy(log(1 - omega) + readProb[i],
                                  log(omega) + readProb[i + nSnp]);
        } else if (refAllele == 1) {
          emitProb[j] = logsum_xy(log(omega) + readProb[i],
                                  log(1 - omega) + readProb[i + nSnp]);
        }
   
        forward[j + (i % 2) * nHap] = 
          logsum_xy(
            forward[j + (i % 2) * nHap],
            vSum[i-1] + tranProbBase[j] + recombProb[snpIdx - 1 + i - 1]) +
          emitProb[j];
      }

      vSum[i] = logsum_array(&forward[(i % 2) * nHap], nHap);
    }
    
    totalVSum[0] = vSum[nSnp - 1];

    free(forward);
    free(emitProb);
    free(vSum);
  
  } else if (nlhs == 4) {
     /* memory efficient forward backward, only keep posterior*/

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nHap, nSnp, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nHap, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nHap, 1, mxREAL);
    
    double* totalVSum = mxGetPr(plhs[0]);
    double* posterior = mxGetPr(plhs[1]);
    double* dist = mxGetPr(plhs[2]);
    double* tranProbBase = mxGetPr(plhs[3]);

    double* backward = malloc_double_array(nHap * 2, 0.0);
    double* previous = malloc_double_array(nHap, 0.0);
    double* previousPlus = malloc_double_array(nHap, 0.0);
    double* emitProb = malloc_double_array(nHap, 0.0);
    double* vSum = malloc_double_array(nSnp, 0.0);
    double previousBase = 0.0;
    
    compute_distance(dist, loc, refloc, hapIdx, nRefHap, nHap);
    compute_transition(tranProbBase, lambda, dist, nHap);

    for (j = 0; j < nHap; j++) {
      int refAllele = getAllele(refhap, 
                                nRefSnp,
                                nRefHap,
                                snpIdx - 1,
                                (int)hapIdx[j] - 1);
      if (refAllele == 0) {
        emitProb[j] = logsum_xy(log(1 - omega) + readProb[0],
                                log(omega) + readProb[0 + nSnp]);
      } else if (refAllele == 1) {
        emitProb[j] = logsum_xy(log(omega) + readProb[0],
                                log(1 - omega) + readProb[0 + nSnp]);
      }
      posterior[j] = tranProbBase[j] + emitProb[j];
    }

    vSum[0] = logsum_array(posterior, nHap);
    for (i = 1; i < nSnp; i++) {
      array_plus_constant(&posterior[(i-1)*nHap], 
                          nHap,
                          nonRecombProb[snpIdx - 1 + i - 1],
                          &posterior[i*nHap]);
    
      for (j = 0; j < nHap; j++) {
        int refAllele = getAllele(refhap, 
                                  nRefSnp,
                                  nRefHap,
                                  snpIdx - 1 + i,
                                  (int)hapIdx[j] - 1); 
        if (refAllele == 0) {
          emitProb[j] = logsum_xy(log(1 - omega) + readProb[i],
                                  log(omega) + readProb[i + nSnp]);
        } else if (refAllele == 1) {
          emitProb[j] = logsum_xy(log(omega) + readProb[i],
                                  log(1 - omega) + readProb[i + nSnp]);
        }
   
        posterior[j + i * nHap] = 
          logsum_xy(
            posterior[j + i * nHap],
            vSum[i-1] + tranProbBase[j] + recombProb[snpIdx - 1 + i - 1]) +
          emitProb[j];
      }

      vSum[i] = logsum_array(&posterior[i * nHap], nHap);
    }
    
    totalVSum[0] = vSum[nSnp - 1];
    
    /* backward */
    for (k = 0; k < nHap; k++) {
      posterior[k + (nSnp - 1) * nHap] -= totalVSum[0];
    }

    for (i = nSnp - 2; i >= 0; i--) {
      for (k = 0; k < nHap; k++) {
        int refAllele = getAllele(refhap, 
                                  nRefSnp,
                                  nRefHap,
                                  snpIdx + i,
                                  (int)hapIdx[k] - 1); 
        if (refAllele == 0) {
          emitProb[k] = logsum_xy(log(1 - omega) + readProb[i + 1],
                                  log(omega) + readProb[i + 1 + nSnp]);
        } else if (refAllele == 1) {
          emitProb[k] = logsum_xy(log(omega) + readProb[i + 1],
                                  log(1 - omega) + readProb[i + 1 + nSnp]);
        }

        previous[k] = backward[k + ((i + 1) % 2) * nHap] +
                      tranProbBase[k] +
                      recombProb[snpIdx - 1 + i] +
                      emitProb[k];

        previousPlus[k] = backward[k + ((i + 1) % 2) * nHap] +
                          nonRecombProb[snpIdx - 1 + i] +
                          emitProb[k];
      }
      previousBase = logsum_array(previous, nHap);
      for (j = 0; j < nHap; j++) {
        backward[j + (i % 2) * nHap] = logsum_xy(previousBase, previousPlus[j]);
        posterior[j + i * nHap] += backward[j + (i % 2) * nHap] - totalVSum[0];
      }
    }

    free(backward);
    free(previous);
    free(previousPlus);
    free(emitProb);
    free(vSum);
  
  } else if (nlhs == 5) {
    /* full computation for everything, but memory consuming */

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nHap, nSnp, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nHap, nSnp, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nHap, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nHap, 1, mxREAL);
    
    double* totalVSum = mxGetPr(plhs[0]);
    double* forward = mxGetPr(plhs[1]);
    double* backward = mxGetPr(plhs[2]);
    double* dist = mxGetPr(plhs[3]);
    double* tranProbBase = mxGetPr(plhs[4]);

    double* previous = malloc_double_array(nHap, 0.0);
    double* previousPlus = malloc_double_array(nHap, 0.0);
    double* emitProb = malloc_double_array(nHap, 0.0);
    double* vSum = malloc_double_array(nSnp, 0.0);
    double previousBase = 0.0;
    
    compute_distance(dist, loc, refloc, hapIdx, nRefHap, nHap);
    compute_transition(tranProbBase, lambda, dist, nHap);

    for (j = 0; j < nHap; j++) {
      int refAllele = getAllele(refhap, 
                                nRefSnp,
                                nRefHap,
                                snpIdx - 1,
                                (int)hapIdx[j] - 1);
      if (refAllele == 0) {
        emitProb[j] = logsum_xy(log(1 - omega) + readProb[0],
                                log(omega) + readProb[0 + nSnp]);
      } else if (refAllele == 1) {
        emitProb[j] = logsum_xy(log(omega) + readProb[0],
                                log(1 - omega) + readProb[0 + nSnp]);
      }
      forward[j] = tranProbBase[j] + emitProb[j];
    }

    vSum[0] = logsum_array(forward, nHap);
    for (i = 1; i < nSnp; i++) {
      array_plus_constant(&forward[(i-1)*nHap], 
                          nHap,
                          nonRecombProb[snpIdx - 1 + i - 1],
                          &forward[i*nHap]);
    
      for (j = 0; j < nHap; j++) {
        int refAllele = getAllele(refhap, 
                                  nRefSnp,
                                  nRefHap,
                                  snpIdx - 1 + i,
                                  (int)hapIdx[j] - 1); 
        if (refAllele == 0) {
          emitProb[j] = logsum_xy(log(1 - omega) + readProb[i],
                                  log(omega) + readProb[i + nSnp]);
        } else if (refAllele == 1) {
          emitProb[j] = logsum_xy(log(omega) + readProb[i],
                                  log(1 - omega) + readProb[i + nSnp]);
        }
   
        forward[j + i * nHap] = 
          logsum_xy(
            forward[j + i * nHap],
            vSum[i-1] + tranProbBase[j] + recombProb[snpIdx - 1 + i - 1]) +
          emitProb[j];
      }

      vSum[i] = logsum_array(&forward[i * nHap], nHap);
    }
    
    totalVSum[0] = vSum[nSnp - 1];

    for (i = nSnp - 2; i >= 0; i--) {
      for (k = 0; k < nHap; k++) {
        int refAllele = getAllele(refhap, 
                                  nRefSnp,
                                  nRefHap,
                                  snpIdx + i,
                                  (int)hapIdx[k] - 1); 
        if (refAllele == 0) {
          emitProb[k] = logsum_xy(log(1 - omega) + readProb[i + 1],
                                  log(omega) + readProb[i + 1 + nSnp]);
        } else if (refAllele == 1) {
          emitProb[k] = logsum_xy(log(omega) + readProb[i + 1],
                                  log(1 - omega) + readProb[i + 1 + nSnp]);
        }

        previous[k] = backward[k + (i + 1) * nHap] +
                      tranProbBase[k] +
                      recombProb[snpIdx - 1 + i] +
                      emitProb[k];

        previousPlus[k] = backward[k + (i + 1) * nHap] +
                          nonRecombProb[snpIdx - 1 + i] +
                          emitProb[k];
      }
      previousBase = logsum_array(previous, nHap);
      for (j = 0; j < nHap; j++) {
        backward[j + i * nHap] = logsum_xy(previousBase, previousPlus[j]);
      }
    }

    free(previous);
    free(previousPlus);
    free(emitProb);
    free(vSum);
  }
}
