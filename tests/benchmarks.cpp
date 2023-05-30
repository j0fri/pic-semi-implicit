#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

int main(int argc, char* argv[]){
    int status = MPI_Init(&argc, &argv);
    if(status != MPI_SUCCESS){
        throw std::runtime_error("MPI initialisation error.");
    }

    //Benchmark 1: landau with Np=10000, scaling Ng
    std::vector<int> Nxs{};
    int Nxi = 16;
    while (Nxi < 10000){
        Nxs.push_back(Nxi);
        Nxi *= 2;
    }
    for(auto Nx: Nxs){
        auto config = preset_configs::landau2D3VX<double>(Nx,10);
        std::stringstream ss;
        ss << "./outputs/benchmark1/semiimplicit/" << Nx << "/";
        config.speciesConfig[0].Np = 10000;
        config.speciesConfig[1].Np = 10000;
        config.saveConfig = {false, false,false, false,false,false,false,
                             false,false,false,true,0.01};
        config.saveConfig.outputFilesDirectory = ss.str();
        config.timeConfig.total = 0.1;
        config.timeConfig.step = 0.01;
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }
    for(auto Nx: Nxs){
        auto config = preset_configs::landau2D3VX<double>(Nx,10);
        std::stringstream ss;
        ss << "./outputs/benchmark1/explicit/" << Nx << "/";
        config.speciesConfig[0].Np = 10000;
        config.speciesConfig[1].Np = 10000;
        config.saveConfig = {false, false,false, false,false,false,false,
                             false,false,false,true,0.01};
        config.saveConfig.outputFilesDirectory = ss.str();
        config.timeConfig.total = 0.1;
        config.timeConfig.step = 0.01;
        config.useExplicitScheme = true;
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }
    //End benchmark 1

    MPI_Finalize();
}