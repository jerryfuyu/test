#ifndef __SOLVELINEARSYSTEM_h_IS_INCLUDED__
#define __SOLVELINEARSYSTEM_h_IS_INCLUDED__
#include <string>
using namespace std;

namespace solveLinearSystem
{
extern int solveLinearSystem(double** A, double* b, double* x, int n, string typeStr, string methodStr, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 1000000);
}
#endif 
