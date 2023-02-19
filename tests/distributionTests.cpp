#include "../models/Distribution.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

template <typename T1>
void outputMatrix(const T1& matrix, std::ostream& output){
	int n = matrix.size();
	int m = matrix[0].size();
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			output << matrix[i][j] << " ";
		}
		output << std::endl;
	}
}

int main(){
    Grid<double,1> grid1{{{{-1,1,100}}}};
    Grid<double,1> grid2{{{{-1,1,100}}}};
    Grid<double,2> grid3{{{{-1,1,100},{-2,2,100}}}};
    Grid<double,1> grid4{{{{-1,1,2}}}};

    std::function f1([](std::array<double,1> input) {return std::exp(-10*input[0]*input[0]);});
    std::function f2([=](std::array<double,1> input) {return (input[0]>0)*f1(input);});
    std::function f3([](std::array<double,2> input){return std::exp(-input[0]*input[0]-input[1]*input[1]);});

    Distribution<double,1> dist1(f1);
    Distribution<double,1> dist2(f2);
    Distribution<double,2> dist3(f3);

    std::ofstream output1("./outputs/distributionTest1");
    std::ofstream output2("./outputs/distributionTest2");
    std::ofstream output3("./outputs/distributionTest3");
    std::ofstream output4("./outputs/distributionTest4");

    outputMatrix(dist1.generate(1000,grid1),output1);
    outputMatrix(dist2.generate(1000,grid2),output2);
    outputMatrix(dist3.generate(1000,grid3),output3);
    outputMatrix(dist1.generate(1000,grid4),output4);
}