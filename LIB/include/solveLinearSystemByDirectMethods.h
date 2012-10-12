#ifndef __SOLVELINEARSYSTEMBYDIRECTMETHODS_h_IS_INCLUDED__
#define __SOLVELINEARSYSTEMBYDIRECTMETHODS_h_IS_INCLUDED__

#include <string>
using namespace std;

namespace solveLinearSystem
{
namespace directMethods
{
//extern int getMethod(string method);

extern void solveLinearSystemByDirectMethods(double** A, double* b, double* x, int n, string method);
}
}
#endif 
