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

//    std::vector<double> dts = {1,0.1,0.01,0.001};
//    int id = 0;
//    for(auto dt: dts){
//        auto config = preset_configs::langmuir(100000,30,5,dt);
//        std::stringstream ss;
//        ss << "./outputs/langmuir/semiimplicit/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//
//    id = 0;
//    for(auto dt: dts){
//        auto config = preset_configs::langmuir(100000,30,5,dt);
//        config.useExplicitScheme = true;
//        std::stringstream ss;
//        ss << "./outputs/langmuir/explicit/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }

//    Single langmuir
    auto config = preset_configs::langmuir(500000,30,3,0.02);
    config.saveConfig.savePositionDistribution = true;
    std::stringstream ss;
    ss << "./outputs/langmuir2/";
    config.saveConfig.outputFilesDirectory = ss.str();
    Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();

    MPI_Finalize();
}