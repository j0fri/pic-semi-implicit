#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

int main(int argc, char* argv[]){    int status = MPI_Init(&argc, &argv);
    if(status != MPI_SUCCESS){
        throw std::runtime_error("MPI initialisation error.");
    }

    std::vector<double> dts = {5,2,1,0.5623,0.316,0.1778,0.1,0.05623,0.0316};
    int id = 0;

    for(auto dt: dts){
        auto config = preset_configs::landau2D3VXCustomBounds<double>(100,3,2.226,10);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.step = dt;

        std::stringstream ss;
        ss << "./outputs/highEnergy/semiimplicit1/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::landau2D3VXCustomBounds<double>(100,3,0.5565,2.5);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.step = dt;

        std::stringstream ss;
        ss << "./outputs/highEnergy/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::landau2D3VXCustomBounds<double>(100,3,0.44296,1.96);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.step = dt;

        std::stringstream ss;
        ss << "./outputs/highEnergy/semiimplicit3/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();


        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;
//
//    for(auto dt: dts){
//        auto config = preset_configs::landau2D3VXCustomBounds<double>(100,3,0.226,1.96);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.step = dt;
//
//        std::stringstream ss;
//        ss << "./outputs/highEnergy/semiimplicit3/" << id << "/";
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
//    for(auto dt: dts){
//        auto config = preset_configs::landau2D3VXCustomBounds<double>(100,3,0.226,1.28);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.step = dt;
//
//        std::stringstream ss;
//        ss << "./outputs/highEnergy/semiimplicit4/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//
//
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//    id = 0;


    MPI_Finalize();
}