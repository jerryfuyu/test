#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "directMethods.h"

using namespace std;

namespace solveLinearSystem
{
namespace directMethods
{
enum Method
{
	Gaussian,
	LU,//for almost all
	Cholesky,
	Thomas,//for tri-diagonal
};

int getMethod(string methodStr)
{
	if(methodStr == "G" || methodStr == "Gaussian" || methodStr == "g" || methodStr == "gaussian")
		return Gaussian;
	if(methodStr == "LU" || methodStr == "lu")
		return LU;
	if(methodStr == "C" || methodStr == "Cholesky" || methodStr == "c" || methodStr == "cholesky")
		return Cholesky;
	if(methodStr == "T" || methodStr == "Thomas" || methodStr == "tri" || methodStr == "triDiag" || methodStr == "thomas" || methodStr == "triDiagonal")
		return Thomas;
	return Gaussian;
}

void solveLinearSystemByDirectMethods(double** A, double* b, double* x, int n, string methodStr)
{
	int method = getMethod(methodStr);
	switch(method)
	{
	case Gaussian:
		solveLinearSystemByGaussianEliminationWithScaledPartialPivoting(A, b, x, n);
		break;
	case Cholesky:
		solveLinearSystemByCholesky(A, b, x, n);
		break;		
	case LU:
		break;	
	case Thomas:
		double* alpha = new double[n];
		double* beta = new double[n];
		double* gamma = new double[n];
		for(int i = 1; i < n; i++)
		{
			beta[i] = A[i][i];
			alpha[i] = A[i][i-1];
			gamma[i-1] = A[i-1][i]; 	
		}
		beta[0] = A[0][0];
		alpha[0] = 0.0;
		gamma[n-1] = 0.0;
		solveLinearSystemByThomasDecom(alpha, beta, gamma, b, x, n);
		delete[] alpha;
		delete[] beta;
		delete[] gamma;
		break;
	}
}
}
}
