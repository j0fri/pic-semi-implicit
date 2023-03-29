#ifndef PIC_SEMI_IMPLICIT_MATH_HELPER_H
#define PIC_SEMI_IMPLICIT_MATH_HELPER_H

#include <vector>

namespace math_helper{
	template <typename T>
	std::vector<std::vector<T>> cartesian(const std::vector<std::vector<T>>& s1, const std::vector<T>& s2);
	
	template <typename T>
	std::vector<std::vector<T>> repeatedCartesian(const std::vector<std::vector<T>>& sets);
	
	template <typename T>
	std::vector<T> linspace(T x1, T x2, unsigned int n);

    //y = A*x, where a is matrix in column-major format
    template <typename T>
    void gemv(unsigned int M, unsigned int N, T alpha, const T* A, unsigned int lda, const T* x,
              unsigned int incx, T* y, unsigned int incy);

}

#include "math_helper.cpp"

#endif //PIC_SEMI_IMPLICIT_MATH_HELPER_H