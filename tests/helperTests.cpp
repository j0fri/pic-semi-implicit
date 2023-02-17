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

int main(){
	std::vector<float> v1 = math_helper::linspace((float)0.0,(float)10.0,21);
	std::vector<float> v2 = math_helper::linspace((float)1.0,(float)10.0,11);
	std::vector<float> v3 = math_helper::linspace((float)1.0,(float)2.0,2);
	
	std::vector<std::vector<float>> vectors{v1, v2, v3};
	
	auto product = math_helper::repeatedCartesian(vectors);
	
	print2DVector(product);
}
