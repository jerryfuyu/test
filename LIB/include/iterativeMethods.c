#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"

namespace solveLinearSystem
{
namespace iterativeMethods
{
int solveByJacobi(double** A, double* b, double* x/*0 is the initial guess*/, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{
// A = D-L-U, then D*x = (L+U)*x + b, so x = inv(D)*(L+U)*x + inv(D)*b	
// so A[i][i] != 0 for each i is a must condition to use Jacobi method
// and if that the method converges globally

	double tmp = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = A[i][i];
		tmp *= tmp;
		if(tmp < abs_eps)
		{
			printf("cannot use Jacobi method to solve the linear system as A[%d][%d] is zero!\n", i, i);
			return -1;
		}
	}
	double *r = new double[n];
	double *y = new double[n];
	for(int i = 0; i < n; i++)
	{
		y[i] = 0.0;
	}
	long iter_count = 0;
	double bnorm = computeVectorMaxNorm(b, n);
	double relErr = 1.0;
	double absErr = bnorm;
	computeMatrixVectorProduct(A, x, y, n);
	while(relErr > rel_eps && absErr > abs_eps)
	{
		iter_count++;
		if(iter_count > max_iter_num)
		{
			printf("iterative number is over %ld!\n", iter_count);
			return -1;
		}
		for(int i = 0; i < n; i++)
		{
			x[i] += (b[i] - y[i])/A[i][i];
		}
		computeMatrixVectorProduct(A, x, y, n);
		computeResidualVector(b, y, r, n);
		absErr = computeResidualNorm(r, n);
		relErr = absErr/bnorm;
		if(iter_count%10000 == 0)
			printf("%.8lf\n", relErr);
	}
	delete [] y;
	delete [] r;
	return iter_count;
}

int solveByGaussSeidel(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{	
// A = D-L-U, then (D-L)*x = U*x + b, so x = inv(D-L)*U*x + inv(D-L)*b	
// so A[i][i] != 0 for each i is a must condition to use Gauss-Seidel method
// and if that the method converges faster than Jacobi, but maybe not converge 
	double tmp = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = A[i][i];
		tmp *= tmp;
		if(tmp < abs_eps)
		{
			printf("cannot use Gauss-Seidel method to solve the linear system as A[%d][%d] is zero!\n", i, i);
			return -1;
		}
	}
	long iter_count = 0;
	double bnorm = computeVectorMaxNorm(b, n);
	double relErr = 1.0;
	double absErr = bnorm;
	double * r = new double[n];
	double sum = 0.0;
	while(relErr > rel_eps && absErr > abs_eps)
	{
		iter_count++;
		if(iter_count > max_iter_num)
		{
			printf("iterative number is over %ld!\n", iter_count);
			return -1;
		}
		for(int i = 0; i < n; i++)
		{	
			sum = b[i];
			for(int j = 0; j < n; j++)
			{
				sum -= A[i][j]*x[j];
			}
			x[i] += sum/A[i][i];
		}
		computeResidualVector(A, b, x, r, n);
		absErr = computeResidualNorm(r, n);
		relErr = absErr/bnorm;
		if(iter_count%10000 == 0)
			printf("%.8lf\n", relErr);
	}
	delete [] r;
	return iter_count;
}


int solveBySOR(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{	
// when we use Gauss-Seildel method : A = D-L-U, then (D-L)*x = U*x + b, so x = inv(D-L)*U*x + inv(D-L)*b to compute the updated x,
// then we use the original x and the updated x to compute the weighted average with the weight omiga which is between 0 and 2,
// then A = 1/omiga*(D-omiga*L) - 1/omiga[(1-omiga)*D + omiga*U] := M - N, and x = inv(M)*N*x + inv(M)*b,
// and if omiga = 1, then the SOR method is the Gauss-Seidel method, 
// the best omiga to A, is 2/(1+sqrt(1-mu*mu)) where mu := rou(A) which means the spectral radius of matrix J, 
// where J means the iterater matrix of the Jacobi method, J := inv(D)*(L+U)
// so A[i][i] != 0 for each i is a must condition to use SOR method
// and if that the method converges faster than GS, but maybe not converge 

/*	double tmp = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = A[i][i];
		tmp *= tmp;
		if(tmp < abs_eps)
		{
			printf("cannot use Gauss-Seidel method to solve the linear system as A[%d][%d] is zero!\n", i, i);
			return -1;
		}
	}
	long iter_count = 0;
	double bnorm = computeVectorMaxNorm(b, n);
	double relErr = 1.0;
	double absErr = bnorm;
	double * r = new double[n];
	double sum = 0.0;
	while(relErr > rel_eps && absErr > abs_eps)
	{
		iter_count++;
		if(iter_count > max_iter_num)
		{
			printf("iterative number is over %ld!\n", iter_count);
			return -1;
		}
		for(int i = 0; i < n; i++)
		{	
			sum = b[i];
			for(int j = 0; j < n; j++)
			{
				sum -= A[i][j]*x[j];
			}
			x[i] += sum/A[i][i];
		}
		computeResidualVector(A, b, x, r, n);
		absErr = computeResidualNorm(r, n);
		relErr = absErr/bnorm;
		if(iter_count%10000 == 0)
			printf("%.8lf\n", relErr);
	}
	delete [] r;
	return iter_count;*/
}

int solveBySSOR(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{	
// we use SOR method two times and symmetric,
// step 1: (D-w*L)*y = [(1-w)*D+w*U]*x + w*b
// step 2: (D-w*U)*x = [(1-w)*D+w*L]*y + w*b
// x = B*x + f, where B := inv(D-w*U)*[(1-w)*D+w*L]*inv(D-w*L)*[(1-w)*D+w*U]
// A = [1/w/(2-w)*(D-w*L)*inv(D)*(D-w*U)]*[I-B] := M - N

/*	double tmp = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = A[i][i];
		tmp *= tmp;
		if(tmp < abs_eps)
		{
			printf("cannot use Gauss-Seidel method to solve the linear system as A[%d][%d] is zero!\n", i, i);
			return -1;
		}
	}
	long iter_count = 0;
	double bnorm = computeVectorMaxNorm(b, n);
	double relErr = 1.0;
	double absErr = bnorm;
	double * r = new double[n];
	double sum = 0.0;
	while(relErr > rel_eps && absErr > abs_eps)
	{
		iter_count++;
		if(iter_count > max_iter_num)
		{
			printf("iterative number is over %ld!\n", iter_count);
			return -1;
		}
		for(int i = 0; i < n; i++)
		{	
			sum = b[i];
			for(int j = 0; j < n; j++)
			{
				sum -= A[i][j]*x[j];
			}
			x[i] += sum/A[i][i];
		}
		computeResidualVector(A, b, x, r, n);
		absErr = computeResidualNorm(r, n);
		relErr = absErr/bnorm;
		if(iter_count%10000 == 0)
			printf("%.8lf\n", relErr);
	}
	delete [] r;
	return iter_count;*/
}

int solveBySpeedestDescent(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{
// it is not a stable method and less useful in practise to solve A*x = b
// the idea is to find the speedest descant direction for the minimization problem (A*x, x) - 2*(b,x)
// the matrix A must be a positive matrix and always be a symmetric matrix

}

int solveByConjugateGradient(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{
//the matrix A must be a positive matrix and always be a symmetric matrix
	long iter_count = 0;
	double bnorm = computeVectorMaxNorm(b, n);
	double absErr = bnorm;
	double relErr = 1.0;
	double *r = new double[n];
	double *p = new double[n];
	double *tmp = new double[n];
	computeResidualVector(A, b, x, r, n);
	for(int i  = 0; i < n; i++)
		p[i] = r[i];

	double alfa = 0.0;
	double beta = 0.0;
	while(relErr > rel_eps && absErr > abs_eps)
	{
		iter_count++;
		if(iter_count > max_iter_num)
		{
			printf("iterative number is over %ld!\n", iter_count);
			return -1;
		}
		computeMatrixVectorProduct(A, p, tmp, n);
		double inner = computeVectorInnerProduct(r, r, n);
		alfa = inner/computeVectorInnerProduct(p, tmp, n);
		for(int i = 0; i < n; i++)
		{
			x[i] += alfa*p[i];
			r[i] -= alfa*tmp[i];
		}
		absErr = computeResidualNorm(r, n);
		relErr = absErr/bnorm;
		beta = computeVectorInnerProduct(r, r, n)/inner;
		for(int i  = 0; i < n; i++)
			p[i] = r[i] + beta*p[i];

		if(iter_count%10000 == 0)
			printf("%.8lf\n", relErr);
	}

	delete [] r;
	delete [] p;
	delete [] tmp;
	return iter_count;
}


inline void computePreconditionerInverseVectorProduct(double **A, double *r, double *s, double *tmp, int n)
{//using the D-L to be the preconditioning matrix
//also we can use the power(D, -0.5) or D-U to be the preconditioning matrix
	for(int i = 0; i < n; i++)
	{
		tmp[i] = r[i];
		for(int j = 0; j < i; j++)
		{
			tmp[i] -= A[i][j]*tmp[j];
		}
		tmp[i] /= A[i][i];
	}
	for(int i = n-1; i >=0; i--)
	{
		s[i] = r[i];
		for(int j = 0; j < n; j++)
		{
			if(j < i)
				s[i] -= A[i][j]*tmp[j];
			if(j > i)
				s[i] -= A[i][j]*s[j];
		}
		s[i] /= A[i][i];
	}
}

int solveByPreconditionedConjugateGradient(double** A, double* b, double* x, int n, double abs_eps = 1.0E-12, double rel_eps = 1.0E-15, long max_iter_num = 10000000)
{
// the preconditioned conjugate gradient method is ofter used in the solution of large linear systems
// in which the matrix is sparse and positive definite
// and the preconditioning matrix P is approximately equal to L in  the Cholesky factorization A = L*Turn(L)
	double bnorm = computeVectorMaxNorm(b, n);
	double absErr = bnorm;
	double relErr = 1.0;

	double *r = new double[n];
	double *s = new double[n];
	double *d = new double[n];
	double alpha = 0.0;
	double beta = 0.0;
	double preSRInner = 0.0;
	double SRInner = 0.0;
	double dAdProduct = 0.0;
 	
	for(int i = 0; i < n; i++)
	{
		r[i] = b[i];
	}
	//////////b will never be used, so i will use it to be a temporary array///////////////
	computePreconditionerInverseVectorProduct(A, r, s, b, n);
	preSRInner = computeVectorInnerProduct(r, s, n);
	for(int i = 0; i < n; i++)
	{
		d[i] = s[i];
	}
	/*
	computePreconditionerInverseVectorProduct(A, r, s, b, n);
	computeMatrixVectorProduct(A, s, b, n);
	dadProduct = computeVectorInnerProduct(s, b, n);
	alpha = preSRInner/dadProduct;
	for(int i = 0; i < n; i++)
	{
		x[i] += alpha*s[i];
		r[i] -= alpha*b[i];
	}*/
	long iter_count = 0;

	while(relErr > rel_eps && absErr > abs_eps)
	{
		iter_count++;
		if(iter_count > max_iter_num)
		{
			printf("iterative number is over %ld!\n", iter_count);
			return -1;
		}
		computeMatrixVectorProduct(A, d, b, n);
		dAdProduct = computeVectorInnerProduct(d, b, n);
		alpha = preSRInner/dAdProduct;
		for(int i = 0; i < n; i++)
		{
			x[i] += alpha*d[i];
			r[i] -= alpha*b[i];
		}
		computePreconditionerInverseVectorProduct(A, r, s, b, n);
		SRInner = computeVectorInnerProduct(r, s, n);
		beta = SRInner/preSRInner;
		preSRInner = SRInner;
		for(int i = 0; i < n; i++)
		{
			d[i] = s[i] + beta*d[i];
		}
		absErr = computeResidualNorm(r, n);
		relErr = absErr/bnorm;
		if(iter_count%10000 == 0)
			printf("%.8lf\n", relErr);
	}
	delete [] r;
	delete [] s;
	delete [] d;
	return iter_count;
}
}
}
