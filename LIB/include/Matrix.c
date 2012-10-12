#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include "Random.h"
void swap(double* v, int n, int i, int j)
{
	assert(i < n && j < n);
	if(i == j)
		return;
	double tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;
}

void interchangeRow(double** A, int m, int n, int i, int j)
{//interchange two rows of matrix A, which size is m*n
	assert(i < m && j < m);
	if(i == j)
		return;
	double tmp = 0.0;
	for(int k = 0; k < n; k++)
	{
		tmp = A[i][k];
		A[i][k] = A[j][k];
		A[j][k] = tmp;
	}
}

void interchangeRow(double** A, int n, int i, int j)
{
	interchangeRow(A, n, n, i, j);
}

int findMaxNormEntryIndex(double* x, int n, int offset = 0)
{
	assert(offset < n);
	double tmp = 0.0;
	double max = 0.0;
	int index = 0;	
	for(int i = offset; i < n; i++)
	{
		tmp = fabs(x[i]);
		//printf("\n%f\n",tmp);
		if(max < tmp)
		{
			index = i; 
			max = tmp;
		}
	}
	return index;
}

double findMaxNormEntry(double* x, int n, int offset = 0)
{
	int index = findMaxNormEntryIndex(x, n, offset);
	return x[index];
}

double getMaxNormEntrySign(double* x, int n)
{
	double sign = (findMaxNormEntry(x, n) > 0) ? 1.0 : -1.0;
	return sign;
}

double computeVectorMaxNorm(double* x, int n)
{
	double max = 0.0;
	double tmp = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = fabs(x[i]);
		max = (max < tmp)? tmp : max;
	}
	return max;
}

double computeResidualNorm(double *r, int n)
{
	return computeVectorMaxNorm(r, n);
}

/*double comupteResidualNormFor2Dim(double** r, int n)
{
	double tmp = 0.0;
	for(int i = 0 ; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			tmp += r[i][j] * r[i][j];
		}
	}
	return sqrt(tmp) / (n-1);
}*/

void computeMatrixVectorProduct(double **a, double *x, double* y, int n)
{
	for(int i = 0; i < n; i++)
	{
		y[i] = 0.0;
		for(int j = 0; j < n; j++)
		{
			y[i] += a[i][j]*x[j];
		}
	}
}

void computeResidualVector(double **a, double *b, double *x, double *r, int n)
{
	computeMatrixVectorProduct(a, x, r, n);
	for(int i = 0; i < n; i++)
		r[i] = b[i] - r[i];
}

void computeResidualVector(double *b, double *y, double *r, int n)
{	
	for(int i = 0; i < n; i++)
		r[i] = b[i] - y[i];
}

double computeVectorInnerProduct(double *x, double *y, int n)
{
	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += x[i]*y[i];
	return sum;
}

double computeVector2Norm(double* x, int n)
{
	double inner = computeVectorInnerProduct(x, x, n);
	return sqrt(inner);
}

void computeMatrixProduct(double **A, double **B, double **AB, int n)
{
	double sum = 0.0;
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			sum = 0.0;
			for(int k = 0; k < n; k++)
			{
				sum += A[i][k]*B[k][j];
			}
			AB[i][j] = sum;
		}
	}
}

double computeMatrixTrace(double** A, int n)
{
	double trace = 0.0;
	for(int i = 0; i < n; i++)
	{
		trace += A[i][i];
	}
	return trace;
}

double computeMatrixFrobeniusNorm(double** A, int n)
{
	double norm = 0.0;
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			norm += A[i][j] * A[i][j];
		}
	}
	return sqrt(norm);
}

void showVector(double *x, int n)
{
	for(int i = 0; i < n; i++)
		printf("%16.10lf\n", x[i]);
}

void showVector(int *x, int n)
{
	for(int i = 0; i < n; i++)
		printf("%d\n", x[i]);
}

void showMatrix(double **A, int m, int n)
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
			printf("%7.4f ", A[i][j]);
		printf("\n");
	}
	printf("\n");
}

void showMatrix(double** A, int n)
{
	if(n > 20)
		showMatrix(A, 20, 20);
	else
		showMatrix(A, n, n);
}

int countNonzero(double **A, int n)
{
	int cnt = 0;
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(fabs(A[i][j]) < 1.0E-8)
				cnt++;
		}
	}
	return cnt;
}

void computeMatrixTransposeMatrixProduct(double** ATA, double** A, int n)
{
	double sum = 0.0;
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			sum = 0.0;
			for(int k = 0; k < n; k++)
			{
				sum += A[k][i]*A[k][j];
			}
			ATA[i][j] = sum;
		}
	}
}

void generateRandomVector(double*x, int n)
{
	for(int i = 0; i < n; i++)
	{
		x[i] = Random(-10.0 , 10.0);
	}
}

void computeMinusOfTwoVectorsWithWeight(double* x, double* y, double* z, int n, double alpha = 1.0)
{
	for(int i = 0; i < n; i++)
	{
		z[i] = x[i] - alpha * y[i];
	}
}

void computeSumOfTwoVectorsWithWeight(double* x, double* y, double* z, int n, double alpha = 1.0)
{
	for(int i = 0; i < n; i++)
	{
		z[i] = x[i] + alpha * y[i];
	}
}

void copyVector(double* x, double* y, int n)
{
	for(int i = 0; i < n; i++)
	{
		y[i] = x[i];
	}
}

double* generateVector(int n)
{
	double* vec = new double[n];
	return vec;
}

void destroyVector(double* vec)
{
	delete[] vec;
}

double** generateMatrix(int m, int n)
{
	double** A = new double*[m];
	for(int i = 0; i < m; i++)
	{
		A[i] = new double[n];
	}
	return A;
}

double** generateMatrix(int n)
{
	return generateMatrix(n, n);
}

void zeroVector(double* v, int n)
{
	for(int i = 0; i < n; i++)
	{
		v[i] = 0.0;
	}
}

void zeroMatrix(double** A, int m, int n)
{
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			A[i][j] = 0.0;
}

void zeroMatrix(double** A, int n)
{
	zeroMatrix(A, n, n);
}

void destroyMatrix(double** A, int n)
{
	for(int i = 0; i < n; i++)
	{
		delete[] A[i];
	}
	delete[] A;
}
