#ifndef __MATRIX_H_IS_INCLUDED__
#define __MATRIX_H_IS_INCLUDED__
extern void swap(double* v, int n, int i, int j);

extern void interchangeRow(double** A, int m, int n, int i, int j);

extern void interchangeRow(double** A, int n, int i, int j);

extern int findMaxNormEntryIndex(double* x, int n, int offset = 0);

extern double findMaxNormEntry(double* x, int n, int offset = 0);

extern double getMaxNormEntrySign(double* x, int n);

extern double computeVectorMaxNorm(double* x, int n);

extern double computeResidualNorm(double *r, int n);

//extern double comupteResidualNormFor2Dim(double** r, int n);

extern void computeMatrixVectorProduct(double **A, double *x, double* y, int n);

extern void computeMatrixProduct(double **A, double **B, double **AB, int n);

extern void computeResidualVector(double **A, double *b, double *x, double *r, int n);

extern void computeResidualVector(double *b, double *y, double *r, int n);

extern double computeVectorInnerProduct(double *x, double *y, int n);

extern double computeVector2Norm(double* x, int n); 

extern double computeMatrixTrace(double** A, int n);

extern double computeMatrixFrobeniusNorm(double** A, int n);

extern void showVector(double *x, int n);

extern void showVector(int *x, int n);

extern void showMatrix(double** A, int m, int n);

extern void showMatrix(double **A, int n);

extern int countNonzero(double **A, int n);

extern void computeMatrixTransposeMatrixProduct(double** ATA, double** A, int n);

extern void generateRandomVector(double*x, int n);

extern void computeMinusOfTwoVectorsWithWeight(double* x, double* y, double* z, int n, double alpha = 1.0);

extern void computeSumOfTwoVectorsWithWeight(double* x, double* y, double* z, int n, double alpha = 1.0);

extern void copyVector(double* x, double* y, int n);

extern double* generateVector(int n);

extern void destroyVector(double* vec);

extern double** generateMatrix(int m, int n);

extern double** generateMatrix(int n);

extern void zeroVector(double* v, int n);

extern void zeroMatrix(double** A, int m, int n);

extern void zeroMatrix(double** A, int n);

extern void destroyMatrix(double** A, int n);


#endif 
