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

    //Benchmark 1: scaling with Np
    std::vector<int> Nps{};
    int Np = 1000;
    while (Np < 10000000){
        Nps.push_back(Np);
        Np *= 2;
    }
    for(auto Npi: Nps){
        auto config = preset_configs::landau2D3VX<double>(30,2);
        std::stringstream ss;
        ss << "./outputs/benchmark2/semiimplicit/" << Npi << "/";
        config.speciesConfig[0].Np = Npi/2;
        config.speciesConfig[1].Np = Npi/2;
        config.saveConfig = {false, false,false, false,false,false,false,
                             false,false,false,true,0.01};
        config.saveConfig.outputFilesDirectory = ss.str();
        config.timeConfig.total = 0.1;
        config.timeConfig.step = 0.01;
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }
    for(auto Npi: Nps){
        auto config = preset_configs::landau2D3VX<double>(30,2);
        std::stringstream ss;
        ss << "./outputs/benchmark2/explicit/" << Npi << "/";
        config.speciesConfig[0].Np = Npi/2;
        config.speciesConfig[1].Np = Npi/2;
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

    MPI_Finalize();
}