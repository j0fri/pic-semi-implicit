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

    std::vector<double> Nxs = {8,12,18,27,40,61,91,137,205};
    int id = 0;
//    for(auto Nx: Nxs){
//        auto config = preset_configs::langmuir(50000,Nx,2,0.01);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.total = 1.0;
//        std::stringstream ss;
//        ss << "./outputs/langmuir3/semiimplicit/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//
//
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//    id = 0;
//
//    for(auto Nx: Nxs){
//        auto config = preset_configs::langmuir(200000,Nx,2,0.01);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.total = 1.0;
//        std::stringstream ss;
//        ss << "./outputs/langmuir3/semiimplicit2/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//
//
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//    id = 0;
//
//    for(auto Nx: Nxs){
//        auto config = preset_configs::langmuir(50000,Nx,2,0.001);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.total = 1.0;
//        std::stringstream ss;
//        ss << "./outputs/langmuir3/semiimplicit3/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//
//
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//    id = 0;

    for(auto Nx: Nxs){
        auto config = preset_configs::langmuir(50000,Nx,2,0.01);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/langmuir3/explicit/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;
    for(auto Nx: Nxs){
        auto config = preset_configs::langmuir(50000,Nx,2,0.001);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/langmuir3/explicit1/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto Nx: Nxs){
        auto config = preset_configs::langmuir(200000,Nx,2,0.001);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/langmuir3/explicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto Nx: Nxs){
        auto config = preset_configs::langmuir(50000,Nx,2,0.0001);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/langmuir3/explicit3/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    MPI_Finalize();
}