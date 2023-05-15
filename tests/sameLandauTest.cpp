#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

int main(int argc, char* argv[]){
    auto config = preset_configs::landauFile<double>();
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