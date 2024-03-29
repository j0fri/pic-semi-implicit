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

    std::vector<double> Factors = {1,1.7783,3.162,5.623,10,17.783,31.623,56.234};
    int id = 0;

//    for(auto Factor: Factors){
//        int Nx = (int)(30*Factor*0.3162);
//        int Np = 10000;
//        double dt = 0.01;
//        auto config = preset_configs::langmuir(Np,Nx,5,dt);
//        config.saveConfig.saveAllTimes = true;
//        config.timeConfig.total = 0.1;
//        std::stringstream ss;
//        ss << "./outputs/hotspots/semiimplicit1/" << id << "/";
//        config.saveConfig.outputFilesDirectory = ss.str();
//        Simulation<double,2,3> sim(config);
//        sim.initialise();
//        sim.run();
//        id += 1;
//    }
//    id = 0;

    for(auto Factor: Factors){
        int Nx = (int)(30*Factor);
        int Np = 100000;
        double dt = 0.01;
        auto config = preset_configs::langmuir(Np,Nx,5,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/hotspots/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto Factor: Factors){
        int Nx = (int)(30*Factor*3.162);
        int Np = 1000000;
        double dt = 0.01;
        auto config = preset_configs::langmuir(Np,Nx,5,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/hotspots/semiimplicit3/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto Factor: Factors){
        int Nx = (int)(30*Factor*10);
        int Np = 10000000;
        double dt = 0.01;
        auto config = preset_configs::langmuir(Np,Nx,5,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/hotspots/semiimplicit4/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;



    MPI_Finalize();
}