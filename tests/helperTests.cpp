#include <iostream>
#include <vector>
#include "../helpers/math_helper.h"

template<typename T>
void printVector(const T& v){
	for(int i = 0; i < (int)v.size()-1; ++i){
		std::cout << v[i] << " ";
	}
	std::cout << v.back() << std::endl;
}

template<typename T>
void print2DVector(const T& matrix){
	int n = matrix.size();
	int m;
	for(int i = 0; i < n; ++i){
		m = matrix[i].size();
		for(int j = 0; j < m; ++j){
			std::cout << matrix[i][j];
			if (j != m - 1){
				std::cout << " ";
			}
		}
		std::cout << std::endl;
	}
}

//Prints M by N matrix stored in column-major format
template<typename T>
void printMatrix(unsigned int M, unsigned int N, unsigned int lda, T* A){
    for(unsigned int i = 0; i < M; ++i){
        for(unsigned int j = 0; j < N; ++j){
            std::cout << A[i+j*lda] << " ";
        }
        std::cout << std::endl;
    }
}

int main(){
	std::vector<float> v1 = math_helper::linspace((float)0.0,(float)10.0,21);
	std::vector<float> v2 = math_helper::linspace((float)1.0,(float)10.0,11);
	std::vector<float> v3 = math_helper::linspace((float)1.0,(float)2.0,2);
	
	std::vector<std::vector<float>> vectors{v1, v2, v3};
	
	auto product = math_helper::repeatedCartesian(vectors);
	
	print2DVector(product);

    std::cout << std::endl << std::endl << std::endl;

    int N = 3;
    auto* A = new double[N*N];
    for(int i = 0; i < N*N; ++i){
        A[i] = i;
    }
    std::cout << "A_obsolete (lda=1): " << std::endl;
    printMatrix(N,N,N,A);

    auto* x = new double[N];
    for(int i = 0; i < N; ++i){
        x[i] = N-i;
    }
    std::cout << "x (incx=1): " << std::endl;
    printMatrix(N,1,N,x);

    auto* y = new double[N];
    math_helper::gemv(N,N,1.0,A,N,x,1,y,1);
    std::cout << "A_obsolete*x: " << std::endl;
    printMatrix(N,1,N,y);

    std::cout << "N-1 by N-1 submatrix of A_obsolete(lda=N) * x(incx=2): " << std::endl;
    auto* y2 = new double[N-1];
    math_helper::gemv(N-1,N-1,1.0,A,N,x,2,y2,1);
    printMatrix(N-1,1,N,y2);

}
