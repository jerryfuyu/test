#ifndef __CONSTANTFUNCTION_H_IS_INCLUDED__
#define __CONSTANTFUNCTION_H_IS_INCLUDED__

namespace constantFunction
{
	typedef double (*ptrFun)(double* x, int n);

	typedef double (*ptrODEFun)(double t, double* y, int n);
//	typedef int (*ptr_i_fun)(double* x, int n);
//	typedef template<typename T> T (*ptr_fun)(T* x, int n);
}

#endif
