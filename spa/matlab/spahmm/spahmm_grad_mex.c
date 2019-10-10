#include "mex.h"
#include "spahmm_common.h"
#include "util.h"

/* *
 * [lamGrad, locGrad] = 
 *    spahmm_fb_mex(loc, readProb, Hapref, snpIdx, hapIdx, Param)
 *
 * This function performs forward-backward algorithm in a stochastic way to
 * compute gradient
 * */

void compute_rawLamGrad(double* lamGrad,
                        const double* dist,
                        const double* hapIdx,
                        const double nHap,
                        const double* tranProbBase) {

  int i;
  double baseGrad = 0.0;
  for (i = 0; i < nHap; i++) {
    baseGrad += exp(tranProbBase[i]) * dist[i];
  }

  for (i = 0; i < nHap; i++) {
    lamGrad[i] = dist[i] - baseGrad;
  }
}

void compute_rawLocGrad(double* locGrad,
                        const double* dist, 
                        const double* loc, 
                        const double* refloc, 
                        const double* hapIdx, 
                        const int nRefHap,
                        const int nHap,
                        const double* tranProbBase) {

  int i;
  for (i = 0; i < nHap; i++) {
    locGrad[2 * i] = (loc[0] - refloc[(int)hapIdx[i] - 1 + 0]) / dist[i];
    locGrad[2 * i + 1] = (loc[1] - refloc[(int)hapIdx[i] - 1 + nRefHap]) / dist[i];
  }

  double baseGrad[2] = {0.0, 0.0};
  for (i = 0; i < nHap; i++) {
    baseGrad[0] += exp(tranProbBase[i]) * locGrad[2 * i];
    baseGrad[1] += exp(tranProbBase[i]) * locGrad[2 * i + 1];
  }
  
  for (i = 0; i < nHap; i++) {
    locGrad[2 * i] -= baseGrad[0];
    locGrad[2 * i + 1] -= baseGrad[1];
  }
}

void sampling(int* randIdx, int nSample, int nHap) {

  int i, j;
  for (i = 0, j = 0; i < nHap && j < nSample; ++i) {
    int ri = nHap - i;
    int rj = nSample - j;
    if (rand() % ri < rj) {
      randIdx[j++] = i;
    }
  }
  
  if (j != nSample) {
    mexErrMsgTxt("Knuth algorithm wrong\n");
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
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

  /* memory efficient forward backward, only keep forward*/

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(2, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  
  double* lamGrad = mxGetPr(plhs[0]);
  double* locGrad = mxGetPr(plhs[1]);
  double* tVsum = mxGetPr(plhs[2]);
  lamGrad[0] = 0.0;
  locGrad[0] = 0.0;
  locGrad[1] = 0.0;

  double* dist = malloc_double_array(nHap, 0.0);
  double* tranProbBase = malloc_double_array(nHap, 0.0);
  double* forward = malloc_double_array(nHap * nSnp, 0.0);
  double* backward = malloc_double_array(nHap * 2, 0.0);
  double* previous = malloc_double_array(nHap, 0.0);
  double* previousPlus = malloc_double_array(nHap, 0.0);
  double* emitProb = malloc_double_array(nHap, 0.0);
  double* vSum = malloc_double_array(nSnp, 0.0);
  double previousBase = 0.0;
  double totalVSum = 0.0;
  int* randIdx = malloc_int_array(nSample, 0);
  double* rawLamGrad = malloc_double_array(nHap, 0.0);
  double* rawLocGrad = malloc_double_array(nHap * 2, 0.0);
  
  
  compute_distance(dist, loc, refloc, hapIdx, nRefHap, nHap);
  compute_transition(tranProbBase, lambda, dist, nHap);
  compute_rawLocGrad(rawLocGrad, dist, loc, refloc, hapIdx, nRefHap, nHap, tranProbBase);
  compute_rawLamGrad(rawLamGrad, dist, hapIdx, nHap, tranProbBase);

  /* forward */
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
  
  totalVSum = vSum[nSnp - 1];
  tVsum[0] = totalVSum;
  
  /* backward */
  for (i = nSnp - 2; i >= 0; i--) {
    
    sampling(randIdx, nSample, nHap);
    double recombRatio = nonRecombProb[i] - recombProb[i];
  
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

      double tranDiff = tranProbBase[k] + recombProb[snpIdx - 1 + i];
      double tranSame = logsum_xy(tranDiff, nonRecombProb[snpIdx - 1 + i]);

      previous[k] = backward[k + ((i + 1) % 2) * nHap] +
                    tranProbBase[k] +
                    recombProb[snpIdx - 1 + i] +
                    emitProb[k];

      previousPlus[k] = backward[k + ((i + 1) % 2) * nHap] +
                        nonRecombProb[snpIdx - 1 + i] +
                        emitProb[k];
      
      /* aggregate for transition from i to i+1 */
      double weightK = 0.0;
      for (j = 0; j < nSample; j++) {
        int jj = randIdx[j];
        double tranProb = (jj == k) ? tranSame : tranDiff;
        double prob = exp(forward[jj + i * nHap] +
                          backward[k + ((i + 1) % 2) * nHap] +
                          tranProb +
                          emitProb[k] -
                          totalVSum);

        double weight = (jj == k) ? exp(recombRatio - tranProbBase[k]) : 0.0;
        weightK += prob / (1.0 + weight);
      }

      locGrad[0] += weightK * rawLocGrad[2 * k];
      locGrad[1] += weightK * rawLocGrad[2 * k + 1];
      lamGrad[0] += weightK * rawLamGrad[k];
    }

    /* compute backward for i */
    previousBase = logsum_array(previous, nHap);
    for (j = 0; j < nHap; j++) {
      backward[j + (i % 2) * nHap] = logsum_xy(previousBase, previousPlus[j]);
    }
  }

  locGrad[0] = - lambda * locGrad[0];
  locGrad[1] = - lambda * locGrad[1];
  lamGrad[0] = - lamGrad[0];

  locGrad[0] = locGrad[0] * nHap / nSample;
  locGrad[1] = locGrad[1] * nHap / nSample;
  lamGrad[0] = lamGrad[0] * nHap / nSample;

  free(dist);
  free(tranProbBase);
  free(forward);
  free(backward);
  free(previous);
  free(previousPlus);
  free(emitProb);
  free(vSum);
  free(randIdx);
  free(rawLamGrad);
  free(rawLocGrad);
}
