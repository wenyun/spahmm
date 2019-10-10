#include <float.h>
#include <math.h>

/* logsum of two scalar */
double logsum_xy(double logX, double logY) {
  if (logX > logY) {
    return logX + log(1.0 + exp(logY - logX));
  } else {
    return logY + log(1.0 + exp(logX - logY));
  }
}

/* logsum all (elements * weight) in an array */
double logsum_array_weight(const double* LogX, const int N, const double weight) {
  int i;
  double LogSum = -DBL_MAX;    
  for (i = 0; i < N; i++) {
    LogSum = logsum_xy(LogSum, LogX[i] * weight);
  }
  return LogSum; 
}

/* logsum all elements in an array */
double logsum_array(const double* LogX, const int N) {
  return logsum_array_weight(LogX, N, 1.0);
}

/* for every element in an array, compute logsum(LogX[i], LogA) */
void logsum_constant(double *LogX, int N, double LogA, double* Output) {
  int i;
  for (i = 0; i < N; i++) {
    Output[i] = logsum_xy(LogX[i], LogA);
  }
}
