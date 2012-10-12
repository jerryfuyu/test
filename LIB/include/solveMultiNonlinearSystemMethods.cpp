#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "solveNonlinearSystemConstants.h"
#include "constantFunction.h"
#include "solveLinearSystem.h"
#include "Matrix.h"

using namespace std;
//using namespace solveLinearSystem;
using namespace constantFunction;

namespace solveNonlinearSystem
{
using namespace nonlinearSystemConstants;

double computeFunValue(ptrFun& pFun, double* x, int n)
{
	return (*pFun)(x, n);
}

void computeFunArrayValue(ptrFun* pFuns, double* x, double* f, int n)
{
	for(int i = 0; i < n; i++)
		f[i] = computeFunValue(pFuns[i], x, n);
}

void computeDerFunMatrixValue(ptrFun** pFunMat, double* x, double** df, int n)
{
	for(int i = 0; i < n; i++)
		computeFunArrayValue(pFunMat[i], x, df[i], n);
}

void computeDerFunValueWithApproximating(ptrFun* pFuns, double* x, double** df, int n, double delta = 1.0E-8)
{
	double* tmp = new double[n];
	for(int i = 0; i < n; i++)
		tmp[i] = x[i];
	for(int j = 0; j < n; j++)
	{
		tmp[j] += delta;
		for(int i = 0; i < n; i++)
		{
			df[i][j] = (computeFunValue(pFuns[i], tmp, n) - computeFunValue(pFuns[i], x, n))/delta;
		}
		tmp[j] -= delta;
	}
	delete[] tmp;
}

int solveNonlinearSystemByNewtonWithDerivatives(ptrFun* pFuns, ptrFun** pDerFunMat, double* x/*the initial guess*/, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, int max_iter_num = 10000000)
{
	double* f = new double[n];
	computeFunArrayValue(pFuns, x, f, n);
	double fNorm = computeVectorMaxNorm(f, n);
	if(fNorm < abs_eps)
		return 0;
	int iterCount = 0;
	double absErr = fNorm;
	double relErr = 1.0;
	double** df = new double* [n];
	for(int i = 0; i < n; i++)
		df[i] = new double[n];
	double* y = new double[n];
	while(relErr > rel_eps && absErr > abs_eps)
	{
		computeDerFunMatrixValue(pDerFunMat, x, df, n);
		solveLinearSystem::solveLinearSystem(df, f, y, n, type, methodStr, abs_eps, rel_eps, max_iter_num); 
		iterCount++;		
		if(iterCount >= max_iter_num)
		{
			printf("iterative number is over %d times!\n", max_iter_num);
			return -1;
		}
		for(int i = 0; i < n; i++)
			x[i] -= y[i];
		//showVector(x, n);
		computeFunArrayValue(pFuns, x, f, n);
		absErr = computeVectorMaxNorm(f, n);
		relErr = absErr/fNorm;
	}
	delete[] y;
	delete[] f;
	for(int i = 0; i < n; i++)
		delete[] df[i];
	delete[] df;	
	return iterCount;
}

int solveNonlinearSystemByNewtonWithApproximatingDerivatives(ptrFun* pFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, double delta = 1.0E-8, int max_iter_num = 10000000)
{
	double* f = new double[n];
	computeFunArrayValue(pFuns, x, f, n);
	double fNorm = computeVectorMaxNorm(f, n);
	if(fNorm < abs_eps)
		return 0;
	int iterCount = 0;
	double absErr = fNorm;
	double relErr = 1.0;
	double** df = new double* [n];
	for(int i = 0; i < n; i++)
		df[i] = new double[n];
	double* y = new double[n];
	while(relErr > rel_eps && absErr > abs_eps)
	{
		computeDerFunValueWithApproximating(pFuns, x, df, n, delta);
		solveLinearSystem::solveLinearSystem(df, f, y, n, type, methodStr, abs_eps, rel_eps, max_iter_num); 
		iterCount++;
		if(iterCount >= max_iter_num)
		{
			printf("iterative numbers is over %d times!\n", max_iter_num);
			return -1;
		}
		for(int i = 0; i < n; i++)
			x[i] -= y[i];
		computeFunArrayValue(pFuns, x, f, n);
		absErr = computeVectorMaxNorm(f, n);
		relErr = absErr/fNorm;
	}
	delete[] y;
	delete[] f;
	for(int i = 0; i < n; i++)
		delete[] df[i];
	delete[] df;
	return iterCount;
}

double getSquareSum(ptrFun* pFuns, double* x, int n)
{
	double sum = 0.0;
	double tmp = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = computeFunValue(pFuns[i], x, n);
		sum += tmp*tmp;
	}
	return sum;
}

void computeGradientWithDerivatives(ptrFun* pFuns, ptrFun** pGradFuns, double* x, double* grad, int n)
{
	double sum = 0.0;
	for(int j = 0; j < n; j++)
	{
		sum = 0.0;
		for(int i = 0; i < n; i++)
		{
			sum += 2.0*computeFunValue(pFuns[i], x, n)*computeFunValue(pGradFuns[j][i], x, n);
		}
		grad[j] = sum;
	}
}

int solveNonlinearSystemBySpeedestDescentWithDerivatives(ptrFun* pFuns, ptrFun** pGradFunMat, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-2, double rel_eps = 1.0E-4, int max_iter_num = 100)
{
// the speedest descent method was presented as a way to obtain good initial approximations.
// the idea is to minimize the sum of f[i]*f[i]
// and if the gradient vector is equal to zero vector, the method can not make sense
	double* grad = new double[n];
	double* minus = new double[n];
	double alpha = 1.0;	
	double g0 = getSquareSum(pFuns, x, n);
	double absErr = g0;	
	double relErr = 1.0;
	int iterCount = 0;
	double gErr = 0.0;
	int alphaCount = 0;
	double grad2Norm = 0.0;
	while(absErr > abs_eps && relErr > rel_eps)
	{
		iterCount++;
		if(iterCount >= max_iter_num)
		{
			printf("iterative numbers is over %d times!\n", max_iter_num);
			return -1;
		}
		//computeGradientWithDerivatives(pFuns, pGradFuns, x, grad, n);		
		//computeMinusOfTwoVectors(x, grad, minus, n, alpha);
		gErr = absErr+1.0;//getSquareSum(pFuns, minus, n);
		alphaCount = 0;
		while(gErr > absErr/2.0 && alphaCount < 10)
		{		
			computeGradientWithDerivatives(pFuns, pGradFunMat, x, grad, n);	
			grad2Norm = computeVector2Norm(grad, n);
			//printf("%f\n",grad2Norm);
			if(grad2Norm < rel_eps)
			{	
				printf("the gradient vector is equal to zero vector, the speedest descent method will not make sense, %f\n", grad2Norm);
				return -1;
			}
			computeMinusOfTwoVectorsWithWeight(x, grad, minus, n, alpha);
			gErr = getSquareSum(pFuns, minus, n);
			alpha /= 2.0;
			alphaCount++;
			//printf("alpha = %.12f\n", alpha);
		}
		absErr = gErr;
		relErr = absErr / g0;	
		//printf("%.16f, %.16f\n", absErr, relErr);
		alpha = 1.0;
		for(int i = 0; i < n; i++)
		{
			x[i] = minus[i];
		}
		//showVector(x, n);
	}
	delete[] grad;
	delete[] minus;
	return iterCount;
}

int solveNonlinearSystemByNewtonAndSpeedestDescentWithDerivatives(ptrFun* pFuns, ptrFun** pDerFunMat, double* x/*the initial guess for speedest descent method*/, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, int max_iter_num = 10000000)
{
	int iterSDM = solveNonlinearSystemBySpeedestDescentWithDerivatives(pFuns, pDerFunMat, x, n);
	int count = 0;
	while(iterSDM == -1 && count++ < 10)
	{
		generateRandomVector(x, n);
		iterSDM = solveNonlinearSystemBySpeedestDescentWithDerivatives(pFuns, pDerFunMat, x, n);
	}
	int iterNewton = solveNonlinearSystemByNewtonWithDerivatives(pFuns, pDerFunMat, x, n);
	printf("iterater number for speedest descent method is %d,\nand for newton method is %d\n", iterSDM, iterNewton);
	return iterNewton;
}

void computeGradientWithApproximatingDerivatives(double* f, double** df, double* grad, int n)
{
	for(int j = 0; j < n; j++)
	{
		for(int i = 0; i < n; i++)
		{
			grad[j] = 2.0*f[i]*df[i][j];
		}
	}
}

int solveNonlinearSystemBySpeedestDescentWithApproximatingDerivatives(ptrFun* pFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-2, double rel_eps = 1.0E-4, double delta = 1.0E-8, int max_iter_num = 100)
{
	double* grad = new double[n];
	double* minus = new double[n];
	double* f = new double[n];
	/*double** df = new double*[n];
	for(int i = 0; i < n; i++)
	{
		df[i] = new double[n];
	}*/
	double** df = generateMatrix(n);
	double alpha = 1.0;
	computeFunArrayValue(pFuns, x, f, n);//f[i]
	double g0 = computeVectorInnerProduct(f, f, n);//sum of f*f
	double absErr = g0;
	double relErr = 1.0;
	int iterCount = 0;
	double gErr = 0.0;
	int alphaCount = 0;
	double grad2Norm = 0.0;
	while(absErr > abs_eps && relErr > rel_eps)
	{
		iterCount++;
		if(iterCount >= max_iter_num)
		{
			printf("iterative numbers is over %d times!\n", max_iter_num);
			return -1;
		}
		//computeGradientWithDerivatives(pFuns, pGradFuns, x, grad, n);		
		//computeMinusOfTwoVectors(x, grad, minus, n, alpha);
		gErr = absErr+1.0;//getSquareSum(pFuns, minus, n);
		alphaCount = 0;
		computeDerFunValueWithApproximating(pFuns, x, df, n, delta);
		computeGradientWithApproximatingDerivatives(f, df, grad, n);
		while(gErr > absErr/2.0 && alphaCount < 10)
		{	
			grad2Norm = computeVector2Norm(grad, n);
			//printf("%f\n",grad2Norm);
			if(grad2Norm < rel_eps)
			{	
				printf("the gradient vector is equal to zero vector, the speedest descent method will not make sense, %f\n", grad2Norm);
				return -1;
			}
			computeMinusOfTwoVectorsWithWeight(x, grad, minus, n, alpha);
			gErr = getSquareSum(pFuns, minus, n);
			alpha /= 2.0;
			alphaCount++;
			//printf("alpha = %.12f\n", alpha);
		}
		absErr = gErr;
		relErr = absErr / g0;	
		//printf("%.16f, %.16f\n", absErr, relErr);
		alpha = 1.0;
		for(int i = 0; i < n; i++)
		{
			x[i] = minus[i];
		}
		computeFunArrayValue(pFuns, x, f, n);
		//showVector(x, n);
	}
	delete[] grad;
	delete[] minus;
	delete[] f;
	/*for(int i = 0; i < n; i++)
	{
		delete[] df[i];
	}
	delete[] df;*/
	destroyMatrix(df, n);
	return iterCount;
}

int solveNonlinearSystemByNewtonAndSpeedestDescentWithApproximatingDerivatives(ptrFun* pFuns, double* x, int n, string type = "D", string methodStr = "G", double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, double delta = 1.0E-8, int max_iter_num = 100)
{
	int iterSDM = solveNonlinearSystemBySpeedestDescentWithApproximatingDerivatives(pFuns, x, n);
	int count = 0;
	while(iterSDM == -1 && count++ < 10)
	{
		generateRandomVector(x, n);
		iterSDM = solveNonlinearSystemBySpeedestDescentWithApproximatingDerivatives(pFuns, x, n);
	}
	showVector(x, n);
	int iterNewton = solveNonlinearSystemByNewtonWithApproximatingDerivatives(pFuns, x, n);
	printf("iterater number for speedest descent method is %d,\nand for newton method is %d\n", iterSDM, iterNewton);
	return iterNewton;
}

}
