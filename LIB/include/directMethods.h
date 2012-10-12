#ifndef __DIRECTMETHODS_h_IS_INCLUDED__
#define __DIRECTMETHODS_h_IS_INCLUDED__
namespace solveLinearSystem
{
namespace directMethods
{ 
extern void solveLinearSystemWithUpperMatrix(double** A, double* b, double* x, int n);

/*extern void solveLinearSystemByGaussianElimination(double** A, double* b, double* x, int n);

extern void solveLinearSystemByGaussianEliminationWithBackwardSubstitution(double** A, double* b, double* x, int n);

extern void solveLinearSystemByGaussianEliminationWithPartialPivoting(double** A, double* b, double* x, int n);*/

extern void solveLinearSystemByGaussianEliminationWithScaledPartialPivoting(double** A, double* b, double* x, int n);

extern void CholeskyDecom(double ** a, int n);

extern void CholeskyY(double** a, double* b, double* y, int n);

extern void CholeskyX(double** a, double* y, double* x, int n);

extern void solveLinearSystemByCholesky(double** A, double* b, double * x, int n);

extern void solveLinearSystemByThomasDecom(double* a, double* b, double* c, double* y, double* x, int n);

extern void computeBetaWithThomasDecom(double* a, double* b, double* c, double* beta, int n);

extern void ThomasDecom(double* a, double* b, double* c, double* beta, double* y, double* x, int n);
}
}
#endif 
