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

    //Benchmark 3: energy loss with dt
    std::vector<int> Nts{};
    int Nt = 10;
    while (Nt < 2000000){
        Nts.push_back(Nt);
        Nt *= 3;
    }
    for(auto Nti: Nts){
        double dt = 10.0/Nti;
        auto config = preset_configs::landau2D3VX<double>(30,2);
        std::stringstream ss;
        ss << "./outputs/benchmark3/semiimplicit/" << Nti << "/";
        config.speciesConfig[0].Np = 10000;
        config.speciesConfig[1].Np = 10000;
        config.saveConfig = {false, false,false, false,true,false,false,
                             true,false,false,true,std::max(dt,0.001)};
        config.saveConfig.outputFilesDirectory = ss.str();
        config.timeConfig.total = 1;
        config.timeConfig.step = dt;
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }
    for(auto Nti: Nts){
        double dt = 10.0/Nti;
        auto config = preset_configs::landau2D3VX<double>(30,2);
        std::stringstream ss;
        ss << "./outputs/benchmark3/explicit/" << Nti << "/";
        config.speciesConfig[0].Np = 10000;
        config.speciesConfig[1].Np = 10000;
        config.saveConfig = {false, false,false, false,true,false,false,
                             true,false,false,true,std::max(dt,0.001)};
        config.saveConfig.outputFilesDirectory = ss.str();
        config.timeConfig.total = 1;
        config.timeConfig.step = dt;
        config.useExplicitScheme = true;
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }

    MPI_Finalize();
}