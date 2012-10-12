#ifndef __SOLVELINEARSYSTEMBYITERATIVEMETHODS_h_IS_INCLUDED__
#define __SOLVELINEARSYSTEMBYITERATIVEMETHODS_h_IS_INCLUDED__

#include <string>
using namespace std;

namespace solveLinearSystem
{
namespace iterativeMethods
{
extern int solveLinearSystemByIterativeMethods(double** A, double* b, double* x, int n, string methodStr, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000);
}
}
#endif 
