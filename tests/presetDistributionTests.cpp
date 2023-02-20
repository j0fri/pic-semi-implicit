#include "../helpers/preset_distributions.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

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

//TODO: test graphically
int main(){
    int Nt = 4;

    Grid<double,1> grid1D{{{{-1,1,100}}}};
    Grid<double,2> grid2D{{{{-1,1,100},{-2,2,100}}}};
    Grid<double,3> grid3D{{{{-1,1,100},{-2,2,100},{-1,1,100}}}};

    //1D boltzmann
    auto dist1 = preset_distributions::Boltzmann<double,1>(10,1,1);
    //3D boltzmann
    auto dist2 = preset_distributions::Boltzmann<double,3>(10,1,1);
    //2D step
    auto dist3 = preset_distributions::Step<double,2>(0,0.5,false);
    //Sine+uniform with step
    auto dist4 = dist3*(preset_distributions::Sin<double,2>(1,1,M_PI,0) + preset_distributions::Uniform<double,2>(1));

    std::vector<std::ofstream> outputs(Nt);
    for(int i = 0; i < Nt; ++i){
        outputs[i] = std::ofstream("./outputs/presetDistributionTest" + std::to_string(i+1));
    }

    outputMatrix(dist1.generate(1000,grid1D),outputs[0]);
    outputMatrix(dist2.generate(1000,grid3D),outputs[1]);
    outputMatrix(dist3.generate(1000,grid2D),outputs[2]);
    outputMatrix(dist4.generate(1000,grid2D),outputs[3]);
}