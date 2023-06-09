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

//    std::vector<double> dts = {0.0316,0.01,0.00316,0.001,0.000316,0.0001};
    std::vector<double> dts = {0.5623,0.316,0.1778,0.1,0.05623,0.0316,0.01,0.00316};
    int id = 0;
    for(auto dt: dts){
        auto config = preset_configs::langmuir(50000,30,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 5.0;
        std::stringstream ss;
        ss << "./outputs/langmuir2/semiimplicit/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }

    id = 0;
    for(auto dt: dts){
        auto config = preset_configs::langmuir(500000,300,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 5.0;
        std::stringstream ss;
        ss << "./outputs/langmuir2/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }

    id = 0;
    for(auto dt: dts){
        auto config = preset_configs::langmuir(50000,30,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.useExplicitScheme = true;
        config.timeConfig.total = 5.0;
        std::stringstream ss;
        ss << "./outputs/langmuir2/explicit/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }

    id = 0;
    for(auto dt: dts){
        auto config = preset_configs::langmuir(500000,300,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.useExplicitScheme = true;
        config.timeConfig.total = 5.0;
        std::stringstream ss;
        ss << "./outputs/langmuir2/explicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }

    MPI_Finalize();
}