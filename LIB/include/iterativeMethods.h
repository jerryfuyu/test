#ifndef __ITERATIVEMETHODS_h_IS_INCLUDED__
#define __ITERATIVEMETHODS_h_IS_INCLUDED__
namespace solveLinearSystem
{
namespace iterativeMethods
{
extern int solveByJacobi(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000);

extern int solveByGaussSeidel(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000);

extern int solveByConjugateGradient(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000);

extern int solveByPreconditionedConjugateGradient(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000);
}
}
#endif 
