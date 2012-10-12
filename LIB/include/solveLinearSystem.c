#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "solveLinearSystemByDirectMethods.h"
#include "solveLinearSystemByIterativeMethods.h"
using namespace std;

//type says using direct or iterative methods to solve the linear system 
//method tells which method is applied to solve the linear system 
namespace solveLinearSystem
{
//using namespace directMethods;
//using namespace iterativeMethods;
enum type
{
	Direct,
	Iterative
};

/*enum method
{
	Jacobi,
	GaussSeidel,
	ConjugateGradient,
	PCG,//preconditioned conjugate Gradient
	FFT,//fast fourier transform
	MultiGrid,
};
*/
int getType(string type)
{
//iterative method is the default method
	int tp = Iterative;
//	type = tolower(type);
	if(type == "d" || type == "D" || type == "direct" || type == "Direct")
		tp = Direct;
	return tp;
}
/*
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
}*/

int solveLinearSystem(double** A, double* b, double* x, int n, string typeStr, string methodStr, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 1000000)
{
	int iterCount = -1;
	int type = getType(typeStr);
//	int method = getMethod(methodStr);
	switch(type)
	{
	case Direct: 
		directMethods::solveLinearSystemByDirectMethods(A, b, x, n, methodStr);
		break;
	case Iterative:
		iterCount = iterativeMethods::solveLinearSystemByIterativeMethods(A, b, x, n, methodStr, abs_eps, rel_eps, max_iter_num);
		break;
	}
	return iterCount;
}
}

