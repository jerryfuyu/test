#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include "Matrix.h"

namespace solveLinearSystem
{
namespace directMethods
{ 

void solveLinearSystemWithUpperMatrix(double** A, double* b, double* x, int n)
{
	//printf("solveLinearSystemWithUpperMatrix!\n");
	for(int i = 0; i < n; i++)
		x[i] = 0.0;
	double sum;
	for(int i = n-1; i >= 0; i--)
	{
		sum = 0.0;
		for(int j = i+1; j < n; j++)
			sum += A[i][j] * x[j];
		x[i] = (b[i] - sum) / A[i][i];
	}
}
/*
void solveLinearSystemByGaussianElimination(double** A, double* b, double* x, int n)
{
}

void solveLinearSystemByGaussianEliminationWithBackwardSubstitution(double** A, double* b, double* x, int n)
{
}

void solveLinearSystemByGaussianEliminationWithPartialPivoting(double** A, double* b, double* x, int n)
{
}*/

void solveLinearSystemByGaussianEliminationWithScaledPartialPivoting(double** A, double* b, double* x, int n)
{
	//printf("solveLinearSystemByGaussianEliminationWithScaledPartialPivoting!\n");
	double* rowMaxV = new double[n];
	for(int i = 0; i < n; i++)
	{
		rowMaxV[i] = findMaxNormEntry(A[i], n);
	}
	double* ratioV = new double[n];
	double tmp = 0.0;
	for(int i = 0; i < n-1; i++)//col
	{
		//showMatrix(A, n);
		//showVector(b, n);
		for(int j = i; j < n; j++)	
		{
			ratioV[j] = fabs(A[j][i] / rowMaxV[j]);
		}
		//compute the pivoting index and swap the rows of matrix A and b
		int pivotingIndex = findMaxNormEntryIndex(ratioV, n, i);
		//printf("\n%d\n",pivotingIndex);
		interchangeRow(A, n, i, pivotingIndex);
		swap(b, n, i, pivotingIndex);		
		swap(rowMaxV, n, i, pivotingIndex);
		for(int k = i+1; k < n; k++)
		{
			tmp = A[k][i] / A[i][i];
			for(int j = i+1; j < n; j++)
			{	
				A[k][j] -= tmp * A[i][j];
			}
			b[k] -= tmp * b[i];
			A[k][i] = 0.0;
		}
	}
	//showMatrix(A, n);
	//showVector(b, n);
	delete[] rowMaxV;
	delete[] ratioV;
	solveLinearSystemWithUpperMatrix(A, b, x, n);
}

bool verifySymmetricMatrix(double** A, int n)
{
	for(int i = 0; i < n; i++)
		for(int j = i+1; j < n; j++)
			if(A[i][j] != A[j][i])
				return false;
	return true;
}

//------------Cholesky Decomposition--------------------///
void CholeskyDecom(double ** A, int n)
{
	for(int j = 0; j < n; j++)
	{	
		for(int k = 0; k < j; k++)
			A[j][j] -= A[j][k]*A[j][k];
		A[j][j] = sqrt(A[j][j]);
		for(int i = j+1; i < n; i++)
		{
			for(int k = 0; k < j; k++)
				A[i][j] -= A[i][k]*A[j][k];
			A[i][j] /= A[j][j];
		}
	}
}

void CholeskyY(double** A, double* b, double* y, int n)
{
	for(int i = 0; i < n; i++)
	{
		double sum = b[i];
		for(int k = 0; k < i; k++)
			sum -= A[i][k]*y[k];
		y[i] = sum/A[i][i];
	}
}

void CholeskyX(double** A, double* y, double* x, int n)
{
	for(int i = n-1; i >= 0; i--)
	{
		double sum = y[i];
		for(int k = i+1; k < n; k++)
			sum -= A[k][i]*x[k];
		x[i] = sum/A[i][i];
	}
}

void solveLinearSystemByCholesky(double** A, double* b, double * x, int n)
{
	//printf("solveLinearSystemByCholesky!\n");
	assert(verifySymmetricMatrix(A, n));
	CholeskyDecom(A, n);
	double* y = new double[n];
	for(int i = 0; i < n; i++)
		y[i] = 0.0;
	CholeskyY(A, b, y, n);
	CholeskyX(A, y, x, n);
	delete[] y;
}

//-------------------------Thomas Algorithm---------------------------------------------//
//-------------------------the best method for solving the linear system----------------//
//-------------------------with the tri-diagonal cofficient matrix----------------------//
void solveLinearSystemByThomasDecom(double* a, double* b, double* c, double* y, double* x, int n)
{//y is the righthandside

	//printf("solveLinearSystemByThomasDecom!\n");
/*	showVector(a, n);
	showVector(b, n);
	showVector(c, n);*/
	double tmp = 0.0;
	c[0] /= b[0];
	y[0] /= b[0];
	//for(int i = 1; i < n-2; i++)
        for(int i = 1; i < n-1; i++)
	{
		tmp = b[i] - a[i] * c[i-1];
		c[i] /= tmp;
		y[i] = (y[i] - a[i] * y[i-1]) / tmp;
                //printf("c[%d]=%f\n",i,c[i]);
	}	
	y[n-1] = (y[n-1] - a[n-1] * y[n-2]) / (b[n-1] - a[n-1] * c[n-2]);
	
	x[n-1] = y[n-1];
	for(int i = n-2; i >= 0; i--)
		x[i] = y[i] - c[i] * x[i+1];
}

void computeBetaWithThomasDecom(double* a, double* b, double* c, double* beta, int n)
{
	double tmp = 0.0;
	beta[0] = c[0] / b[0];
        for(int i = 1; i < n-1; i++)
	{
		tmp = b[i] - a[i] * beta[i-1];
		beta[i] = c[i] / tmp;
	}
}

void ThomasDecom(double* a, double* b, double* c, double* beta, double* y, double* x, int n)
{
	y[0] /= b[0];
        for(int i = 1; i < n-1; i++)
	{
		y[i] = (y[i] - a[i] * y[i-1]) * beta[i] / c[i];
	}	
	y[n-1] = (y[n-1] - a[n-1] * y[n-2]) / (b[n-1] - a[n-1] * beta[n-2]);
	
	x[n-1] = y[n-1];
	for(int i = n-2; i >= 0; i--)
		x[i] = y[i] - beta[i] * x[i+1];
}
}
}

