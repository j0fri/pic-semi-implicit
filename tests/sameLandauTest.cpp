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

    auto config = preset_configs::landauFile<double>();
    config.saveConfig.saveElectrostaticPotential = true;
    Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();

    std::vector<std::string> fileNames = {
            "speciesPositionDistribution",
            "electricField",
            "magneticField",
            "fieldEnergy",
            "speciesEnergy",
            "electrostaticPotential"
    };

    int processId, numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    if(processId == 0){
        for(auto fileName: fileNames){
            std::ifstream refFile("referenceOutputs/landau/" + fileName + ".txt", std::ios::in);
            std::ifstream file("outputs/" + fileName + ".txt", std::ios::in);
            if(output_helper::testSameFileContent<double>(file,refFile, std::pow(10,-12))){
                std::cout << fileName << " files SAME" << std::endl;
            }else{
                std::cout << fileName << " files NOT SAME" << std::endl;
            }
        }
    }

    MPI_Finalize();
}