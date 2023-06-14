#include <iostream>
#include <mpi/mpi.h>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"


double getRuntimePerStep(double dt, int Np, int Ng, bool exp = false){
    auto config = preset_configs::langmuir(Np,Ng/2,2,dt);
    config.saveConfig.outputFilesDirectory = "./outputs/parallelScaling/dummy/";
    config.saveConfig.saveAllTimes = true;
    config.timeConfig.total = dt*10;
    config.useExplicitScheme = exp;

    Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();

    int processId, numProcesses;
    int rankStatus = MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    int sizeStatus = MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    if(rankStatus != MPI_SUCCESS || sizeStatus != MPI_SUCCESS){
        throw std::runtime_error("Could not obtain MPI rank.");
    }

    if(processId > 0){
        return 1;
    }

    std::ifstream dummyRuntimes("./outputs/parallelScaling/dummy/runtime.txt");
    if(!dummyRuntimes.is_open()){
        throw std::runtime_error("Dummy runtime file not open.");
    }

    double initTime = 0;
    double runTime = 0;
    dummyRuntimes >> initTime;
    dummyRuntimes >> runTime;
    dummyRuntimes.close();

    return runTime/10;
}


int main(int argc, char* argv[]){
    int status = MPI_Init(&argc, &argv);
    if(status != MPI_SUCCESS){
        throw std::runtime_error("MPI initialisation error.");
    }

    int processId, numProcesses;
    int rankStatus = MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    int sizeStatus = MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    if(rankStatus != MPI_SUCCESS || sizeStatus != MPI_SUCCESS){
        throw std::runtime_error("Could not obtain MPI rank.");
    }


    //Particle-heavy test:
    {
        double runtimePerStep = getRuntimePerStep(0.1,1000000,90,false);
        int steps = 1.0/runtimePerStep;
        MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);

        auto config = preset_configs::langmuir(1000000,30,3,0.1);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = steps*0.1;
        std::stringstream ss;
        ss << "./outputs/parallelScaling/particleHeavy/" << numProcesses << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }

    //Balanced test:
    {
        double runtimePerStep = getRuntimePerStep(0.1,1000000,3000,false);
        int steps = 1.0/runtimePerStep;
        MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);

        auto config = preset_configs::langmuir(1000000,100,30,0.1);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = steps*0.1;
        std::stringstream ss;
        ss << "./outputs/parallelScaling/balanced/" << numProcesses << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }

    //Grid-heavy test:
    {
        double runtimePerStep = getRuntimePerStep(0.1,25000,5000,false);
        int steps = 5.0/runtimePerStep;
        MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);

        auto config = preset_configs::langmuir(25000,100,50,0.1);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = steps*0.1;
        std::stringstream ss;
        ss << "./outputs/parallelScaling/gridHeavy/" << numProcesses << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }

    //Explicit:
    {
        double runtimePerStep = getRuntimePerStep(0.1,1000000,90,false);
        int steps = 5.0/runtimePerStep;
        MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);

        auto config = preset_configs::langmuir(100000,30,3,0.1);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = steps*0.1;
        std::stringstream ss;
        ss << "./outputs/parallelScaling/explicit/" << numProcesses << "/";
        config.saveConfig.outputFilesDirectory = ss.str();

        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();
    }


    MPI_Finalize();
}