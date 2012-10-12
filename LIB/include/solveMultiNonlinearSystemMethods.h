#ifndef __SOLVEMULTINONLINEARSYSTEMMETHODS_H_IS_INCLUDED__
#define __SOLVEMULTINONLINEARSYSTEMMETHODS_H_IS_INCLUDED__
#include <string>
#include "constantFunction.h"
using namespace std;
using namespace constantFunction;

namespace solveNonlinearSystem
{
extern int solveNonlinearSystemByNewtonWithDerivatives(ptrFun* pFuns, ptrFun** pDerFuns, double* x/*the initial guess*/, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, int max_iter_num = 10000000);

extern int solveNonlinearSystemByNewtonWithApproximatingDerivatives(ptrFun* pFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, double delta = 1.0E-8, int max_iter_num = 10000000);

extern int solveNonlinearSystemBySpeedestDescentWithDerivatives(ptrFun* pFuns, ptrFun** pGradFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-3, double rel_eps = 1.0E-6, int max_iter_num = 100);

extern int solveNonlinearSystemByNewtonAndSpeedestDescentWithDerivatives(ptrFun* pFuns, ptrFun** pDerFunMat, double* x/*the initial guess for speedest descent method*/, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, int max_iter_num = 10000000);

extern int solveNonlinearSystemBySpeedestDescentWithApproximatingDerivatives(ptrFun* pFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-2, double rel_eps = 1.0E-4, double delta = 1.0E-8, int max_iter_num = 100);

extern int solveNonlinearSystemByNewtonAndSpeedestDescentWithApproximatingDerivatives(ptrFun* pFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, double delta = 1.0E-8, int max_iter_num = 100);

}
#endif 
