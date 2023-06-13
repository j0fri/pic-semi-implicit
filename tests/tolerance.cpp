#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

int main(int argc, char* argv[]){    int status = MPI_Init(&argc, &argv);
    if(status != MPI_SUCCESS){
        throw std::runtime_error("MPI initialisation error.");
    }

    std::vector<double> dts = {0.5623,0.316,0.1778,0.1,0.05623,0.0316,0.01,0.00316,0.001};
    int id = 0;
    for(auto dt: dts){
        auto config = preset_configs::langmuir(100000,30,3,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        std::stringstream ss;
        ss << "./outputs/tolerance/semiimplicit/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        config.fieldConfig.solverTolerance = 1;
        config.saveConfig.saveSolverSteps = true;

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::langmuir(100000,30,3,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        std::stringstream ss;
        ss << "./outputs/tolerance/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        config.fieldConfig.solverTolerance = 1e-2;
        config.saveConfig.saveSolverSteps = true;

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::langmuir(100000,30,3,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        std::stringstream ss;
        ss << "./outputs/tolerance/semiimplicit3/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        config.fieldConfig.solverTolerance = 1e-5;
        config.saveConfig.saveSolverSteps = true;

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::langmuir(100000,30,3,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        std::stringstream ss;
        ss << "./outputs/tolerance/semiimplicit4/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        config.fieldConfig.solverTolerance = 1e-8;
        config.saveConfig.saveSolverSteps = true;

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

//    for(auto dt: dts){
//        auto config = preset_configs::langmuir(100000,300,3,dt);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.total = 1.0;
//        std::stringstream ss;
//        ss << "./outputs/tolerance/semiimplicit5/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//
//        config.fieldConfig.solverTolerance = 1e-50;
//        config.saveConfig.saveSolverSteps = true;
//
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::langmuir(100000,30,3,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/tolerance/explicit/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }




    MPI_Finalize();
}