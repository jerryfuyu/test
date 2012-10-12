#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "iterativeMethods.h"

using namespace std;

namespace solveLinearSystem
{
namespace iterativeMethods
{
enum Method
{
	Jacobi,
	GaussSeidel,
	ConjugateGradient,
	PCG,//preconditioned conjugate Gradient
	FFT,//fast fourier transform
	MultiGrid,
};

int getMethod(string method)
{
//jacobi method is the default method
	if(method == "J" || method == "Jacobi" || method == "j" || method == "jacobi")
		return Jacobi;
	if(method == "G" || method == "GS" || method == "GaussSeidel" || method == "g" || method == "gs" || method == "gaussSeidel")
		return GaussSeidel;
	if(method == "C" || method == "CG" || method == "ConjugateGradient" || method == "c" || method == "cg" ||  method == "conjugateGradient")
		return ConjugateGradient;
	if(method == "PCG" || method == "pcg" || method == "preconditioned conjugate Gradient")
		return PCG;
	if(method == "FFT" || method == "fft" || method == "fast fourier transform")
		return FFT;
	if(method == "MultiGrid" || method == "M" || method == "MG" || method == "m" || method == "mg" || method == "multigrid" || method == "MultiGrid")
		return MultiGrid;
	return Jacobi;
}

int solveLinearSystemByIterativeMethods(double** A, double* b, double* x, int n, string methodStr, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{
	int method = getMethod(methodStr);
	int iterCount = 0;
	switch(method)
	{
	case Jacobi:
		iterCount = solveByJacobi(A, b, x, n);
		break;
	case GaussSeidel:
		iterCount = solveByGaussSeidel(A, b, x, n);
		break;
	case ConjugateGradient:
		iterCount = solveByConjugateGradient(A, b, x, n);
		break;
	case PCG:
		iterCount = solveByPreconditionedConjugateGradient(A, b, x, n);
		break;
	case FFT:
		break;
	case MultiGrid:
		break;
	}
	return iterCount;	
}
}
}

