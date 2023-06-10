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

    std::vector<double> Powers = {-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5};
    int id = 0;

    double runtimePerStepObj = 1e-3;
    for(auto power: Powers){
        int Np;
        if(power < -3.5){
            Np = 1000;
        }else{
            Np = 100000;
        }
        int Ng = std::max((int)std::sqrt(Np/std::pow(10,power)),2);
        double dt = 0.01;
        auto dummyConfig = preset_configs::langmuir(Np,Ng/2,2,dt);
        dummyConfig.saveConfig.saveAllTimes = true;
        dummyConfig.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio/semiimplicit1/dummy" << id << "/";
        dummyConfig.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(dummyConfig);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "runtime.txt";
        std::ifstream dummyRuntimes(ss2.str());

        if(!dummyRuntimes.is_open()){
            throw std::runtime_error("Dummy runtime file not open.");
        }

        double initTime = 0;
        double runTime = 0;
        dummyRuntimes >> initTime;
        dummyRuntimes >> runTime;
        dummyRuntimes.close();

        double runtimePerStep = runTime/10;
        int actualNp = (int)(Np*runtimePerStepObj/runtimePerStep);
        int actualNg = (int)std::sqrt(actualNp/std::pow(10,power));

        if(actualNp < 10){
            continue;
        }

        auto config = preset_configs::langmuir(actualNp,actualNg/2,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss3;
        ss3 << "./outputs/optratio/semiimplicit1/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss3.str();
        Simulation<double,2,3> sim2(config);
        sim2.initialise();
        sim2.run();

        std::stringstream ss4;
        ss4 << ss3.str();
        ss4 << "Np.txt";
        std::ofstream npfile(ss4.str());
        npfile << actualNp;
        npfile.close();

        id += 1;
    }
    id = 0;

    runtimePerStepObj = 1e-4;
    for(auto power: Powers){
        if(power > 3){
            continue;
        }
        int Np = 10000;
        int Ng = std::max((int)std::sqrt(Np/std::pow(10,power)),2);
        double dt = 0.01;
        auto dummyConfig = preset_configs::langmuir(Np,Ng/2,2,dt);
        dummyConfig.saveConfig.saveAllTimes = true;
        dummyConfig.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio/semiimplicit2/dummy" << id << "/";
        dummyConfig.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(dummyConfig);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "runtime.txt";
        std::ifstream dummyRuntimes(ss2.str());

        if(!dummyRuntimes.is_open()){
            throw std::runtime_error("Dummy runtime file not open.");
        }

        double initTime = 0;
        double runTime = 0;
        dummyRuntimes >> initTime;
        dummyRuntimes >> runTime;
        dummyRuntimes.close();

        double runtimePerStep = runTime/10;
        int actualNp = (int)(Np*runtimePerStepObj/runtimePerStep);
        int actualNg = (int)std::sqrt(actualNp/std::pow(10,power));

        if(actualNp < 10){
            continue;
        }

        auto config = preset_configs::langmuir(actualNp,actualNg/2,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss3;
        ss3 << "./outputs/optratio/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss3.str();
        Simulation<double,2,3> sim2(config);
        sim2.initialise();
        sim2.run();

        std::stringstream ss4;
        ss4 << ss3.str();
        ss4 << "Np.txt";
        std::ofstream npfile(ss4.str());
        npfile << actualNp;
        npfile.close();

        id += 1;
    }
    id = 0;

    runtimePerStepObj = 1e-2;
    for(auto power: Powers){
        int Np = 100000;
        int Ng = std::max((int)std::sqrt(Np/std::pow(10,power)),2);
        double dt = 0.01;
        auto dummyConfig = preset_configs::langmuir(Np,Ng/2,2,dt);
        dummyConfig.saveConfig.saveAllTimes = true;
        dummyConfig.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio/semiimplicit3/dummy" << id << "/";
        dummyConfig.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(dummyConfig);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "runtime.txt";
        std::ifstream dummyRuntimes(ss2.str());

        if(!dummyRuntimes.is_open()){
            throw std::runtime_error("Dummy runtime file not open.");
        }

        double initTime = 0;
        double runTime = 0;
        dummyRuntimes >> initTime;
        dummyRuntimes >> runTime;
        dummyRuntimes.close();

        double runtimePerStep = runTime/10;
        int actualNp = (int)(Np*runtimePerStepObj/runtimePerStep);
        int actualNg = (int)std::sqrt(actualNp/std::pow(10,power));

        if(actualNp < 10){
            continue;
        }

        auto config = preset_configs::langmuir(actualNp,actualNg/2,2,dt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss3;
        ss3 << "./outputs/optratio/semiimplicit3/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss3.str();
        Simulation<double,2,3> sim2(config);
        sim2.initialise();
        sim2.run();

        std::stringstream ss4;
        ss4 << ss3.str();
        ss4 << "Np.txt";
        std::ofstream npfile(ss4.str());
        npfile << actualNp;
        npfile.close();

        id += 1;
    }
    id = 0;


    MPI_Finalize();
}