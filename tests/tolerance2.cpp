#include <iostream>
#include <limits>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

int main(int argc, char* argv[]){    int status = MPI_Init(&argc, &argv);
    if(status != MPI_SUCCESS){
        throw std::runtime_error("MPI initialisation error.");
    }

    std::vector<double> dts = {0.05623,0.0316,0.017778,0.01,0.005623,0.00316};
    int id = 0;

    for(auto dt: dts){
        auto config = preset_configs::langmuir<double>(200000,30,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        std::stringstream ss;
        ss << "./outputs/tolerance2/semiimplicit1/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        config.fieldConfig.solverTolerance = 1e-8;
        config.saveConfig.saveSolverSteps = true;

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    for(auto dt: dts){
        auto config = preset_configs::langmuir<float>(200000,30,2,(float)dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1.0;
        std::stringstream ss;
        ss << "./outputs/tolerance2/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        config.fieldConfig.solverTolerance = 1e-8;
        config.saveConfig.saveSolverSteps = true;

        Simulation<float,2,3> sim(config);
        sim.initialise();
        sim.run();
        id += 1;
    }
    id = 0;

    std::cout << "float limit: " << std::numeric_limits<float>::min() << std::endl;
    std::cout << "double limit: " << std::numeric_limits<double>::min() << std::endl;

    MPI_Finalize();
}