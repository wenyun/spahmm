#include <float.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

/* structure for Array of two dimensions */
typedef struct {
  double* Data;
  int M;
  int N;
} MatrixType;

/* structure for Array of three dimensions */
typedef struct {
  double* Data;
  int M;
  int N;
  int L;
} Matrix3Type;

/* multi nomial distribution random number */
int mnrnd(double* Prob, int N) {

  int i;

  int Index = -1;
  while (Index < 0) {
    double Rand = rand() / (double)RAND_MAX;
    double Sum = 0.0;
    for (i = 0; i < N; i++) {
      Sum += Prob[i];
      if (Sum >= Rand) {
        Index = i + 1;
        break;
      }
    }
  }

  return Index;
}

/* malloc an array of size N with all element set to Value */
double* malloc_double_array(int N, double Value) {
  int i;
  double* Array = Malloc(double, N);
  for (i = 0; i < N; i++) {
    Array[i] = Value;
  }
  return Array;
}

/* malloc an array of size N with all element set to Integer Value */
int* malloc_int_array(int N, int Value) {
  int i;
  int* Array = Malloc(int, N);
  for (i = 0; i < N; i++) {
    Array[i] = Value;
  }
  return Array;
}

/* initialize an existing array by random integer */
void initialize_int_random(int *A, int N) {
  int i;
  for (i = 0; i < N; i++) {
    double Rand = rand() / (double) RAND_MAX;
    if (Rand > 0.5) {
      A[i] = 1;
    } else {
      A[i] = -1;
    }
  }
}

/* initialize an existing array by some integer number */
void initialize_int(int *A, int N, int Value) {
  int i;
  for (i = 0; i < N; i++) {
    A[i] = Value;
  }
}

/* array plus a constant for every element, to another array */
void array_plus_constant(double* Array, int N, double Value, double* Dest) {
  int i;
  for (i = 0; i < N; i++) {
    Dest[i] = Array[i] + Value;
  }
}

/* weighted array plus weighted array, to another array*/
void array_plus_array(double* Array1,
                      double Weight1,
                      double* Array2,
                      double Weight2,
                      double* Dest,
                      int N) {
  int i;
  for (i = 0; i < N; i++) {
    Dest[i] = Array1[i] * Weight1 + Array2[i] * Weight2;
  }
}

/* array copying */
void array_copy(double* Source, int N, double* Dest) {
  int i;
  for (i = 0; i < N; i++) {
    Dest[i] = Source[i];
  }
}

/* normalize the array, to sum to one */
void normalize(double* Prob, int N) {
  int i;
  double Sum = 0;
  for (i = 0; i < N; i++) {
    Sum += Prob[i];
  }
  for (i = 0; i < N; i++) {
    Prob[i] /= Sum;
  }
}

/* get array value, two dimensions */
double get_value(double* Array, int i, int j, int M, int N) {
  if (i >= M || j >= N) {
    printf("Index out of bound\n");
    exit(1);
  }
  return Array[i + j * M];
}

/* get element i, j  */
double get_matrix(MatrixType* Matrix, int i, int j) {
  return Matrix->Data[i + j * Matrix->M];
}

/* set element i, j to be of Value */
void set_matrix(MatrixType* Matrix, int i, int j, double Value) {
  Matrix->Data[i + j * Matrix->M] = Value;
}

/* add element i, j by Value */
void add_matrix(MatrixType* Matrix, int i, int j, double Value) {
  Matrix->Data[i + j * Matrix->M] += Value;
}

/* initialize the matrix to be of Value */
void initialize_matrix(MatrixType* Matrix, double Value) {
  int i;
  int j;
  for (i = 0; i < Matrix->M; i++) {
    for (j = 0; j < Matrix->N; j++) {
      set_matrix(Matrix, i, j, Value);
    }
  }
}

/* malloc matrix struct */
MatrixType* malloc_matrix(int M, int N, double Value) {
  MatrixType* Matrix = Malloc(MatrixType, 1);
  Matrix->Data = Malloc(double, M * N);
  Matrix->M = M;
  Matrix->N = N;
  initialize_matrix(Matrix, Value);
  return Matrix;
}

/* create matrix struct from pointers and dimensions  */
MatrixType* create_matrix(double* Data, int M, int N) {
  MatrixType* Matrix = Malloc(MatrixType, 1);
  Matrix->Data = Data;
  Matrix->M = M;
  Matrix->N = N;
  return Matrix;
}

/* copy matrix data to Dest */
void copy_matrix_data(MatrixType* Matrix, double* Dest) {
  int i;
  int j;
  for (j = 0; j < Matrix->N; j++) {
    for (i = 0; i < Matrix->M; i++) {
      Dest[i + j * Matrix->M] = get_matrix(Matrix, i, j);
    }
  }
}

/* free the Matrix struct only  */
void shallow_free_matrix(MatrixType* Matrix) {
  free(Matrix);
}

/* free the Matrix struct and data both */
void deep_free_matrix(MatrixType* Matrix) {
  free(Matrix->Data);
  free(Matrix);
}

/* get element i, j, k */
double get_matrix3(Matrix3Type* Matrix, int i, int j, int k) {
  return Matrix->Data[i + j * Matrix->M + k * Matrix->M * Matrix->N];
}

/* set element i, j, k to be of Value */
void set_matrix3(Matrix3Type* Matrix, int i, int j, int k, double Value) {
  Matrix->Data[i + j * Matrix->M + k * Matrix->M * Matrix->N] = Value;
}

/* add element i, j, k by Value */
void add_matrix3(Matrix3Type* Matrix, int i, int j, int k, double Value) {
  Matrix->Data[i + j * Matrix->M + k * Matrix->M * Matrix->N] += Value;
}

/* initialize the matrix3 to be of Value */
void initialize_matrix3(Matrix3Type* Matrix, double Value) {
  int i;
  int j;
  int k;
  for (i = 0; i < Matrix->M; i++) {
    for (j = 0; j < Matrix->N; j++) {
      for (k = 0; k < Matrix->L; k++) {
        set_matrix3(Matrix, i, j, k, Value);
      }
    }
  }
}

/* malloc matrix3 struct */
Matrix3Type* malloc_matrix3(int M, int N, int L, double Value) {
  Matrix3Type* Matrix = Malloc(Matrix3Type, 1);
  Matrix->Data = Malloc(double, M * N * L);
  Matrix->M = M;
  Matrix->N = N;
  Matrix->L = L;
  initialize_matrix3(Matrix, Value);
  return Matrix;
}

/* create matrix3 struct from pointers and dimensions  */
Matrix3Type* create_matrix3(double* Data, int M, int N, int L) {
  Matrix3Type* Matrix = Malloc(Matrix3Type, 1);
  Matrix->Data = Data;
  Matrix->M = M;
  Matrix->N = N;
  Matrix->L = L;
  return Matrix;
}

/* free the Matrix struct only  */
void shallow_free_matrix3(Matrix3Type* Matrix) {
  free(Matrix);
}

/* free the Matrix struct and data both */
void deep_free_matrix3(Matrix3Type* Matrix) {
  free(Matrix->Data);
  free(Matrix);
}
