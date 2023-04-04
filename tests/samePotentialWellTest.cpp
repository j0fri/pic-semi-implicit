#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

int main(int argc, char* argv[]){
    auto config = preset_configs::constPotentialWellFile<double>();
    Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();

    std::ifstream posRefFile("referenceOutputs/constPotentialWell/speciesPosition.txt", std::ios::in);
    std::ifstream velRefFile("referenceOutputs/constPotentialWell/speciesVelocity.txt", std::ios::in);
    std::ifstream posFile("outputs/speciesPosition.txt", std::ios::in);
    std::ifstream velFile("outputs/speciesVelocity.txt", std::ios::in);

    if(output_helper::testSameFileContent<double>(posFile,posRefFile)){
        std::cout << "Position files SAME" << std::endl;
    }else{
        std::cout << "Position files NOT SAME" << std::endl;
    }
    if(output_helper::testSameFileContent<double>(velFile,velRefFile)){
        std::cout << "Velocity files SAME" << std::endl;
    }else{
        std::cout << "Velocity files NOT SAME" << std::endl;
    }
}