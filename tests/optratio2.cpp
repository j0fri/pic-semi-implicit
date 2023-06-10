#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

double getRuntimePerStep(double dt, int Np, int Ng){
    std::cout << "Np=" << Np << std::endl;
    auto config = preset_configs::langmuir(Np,Ng/2,2,dt);
    config.saveConfig.outputFilesDirectory = "./outputs/optratio2/dummy/";
    config.saveConfig.saveAllTimes = true;
    config.timeConfig.total = dt*10;

    Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();

    std::ifstream dummyRuntimes("./outputs/optratio2/dummy/runtime.txt");
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

int searchNp(double NpNt, double NpNg2, double totalRuntime){
//    std::vector<double> Nps;
//    std::vector<double> Nps;
//    for(int i = 0; i < 7; ++i){
//        Np
//        Nps.push_back((int)(2*std::pow(10,(double)i)));
//
//        runtimes.push_back(getRuntimePerStep())
//    }
//
//    tk::spline s(X,Y);
//    double x=1.5;
//    double y=s(x);


    int Np = std::max((int)NpNt+1,1000);

    int Np1 = Np;
    int Ng1 = std::max((int)std::sqrt(Np1/NpNg2),2);
    int Nt1 = (double)Np1/NpNt;
    double runtime1 = Nt1*getRuntimePerStep((double)0.1/Nt1, Np1,Ng1);

    int Np2 = std::min(std::max((int)(Np1* totalRuntime/runtime1),100),500000);
    int Ng2 = std::max((int)std::sqrt(Np2/NpNg2),2);
    int Nt2 = (double)Np2/NpNt;
    double runtime2 = Nt2*getRuntimePerStep((double)0.1/Nt2, Np2,Ng2);

    int Np3 = Np2/2;
    int Ng3 = std::max((int)std::sqrt(Np3/NpNg2),2);
    int Nt3 = (double)Np3/NpNt;
    double runtime3 = Nt3*getRuntimePerStep((double)0.1/Nt3, Np3,Ng3);

    Np1 = Np2;
    runtime1 = runtime2;
    Np2 = Np3;
    runtime2 = runtime3;

    for(int i = 0; i < 10; ++i){
        int temp = Np2;
        Np2 = Np1 + (totalRuntime-runtime1)*(Np2-Np1)/(runtime2-runtime1);
        Np2 = std::max(Np2,2);
        Np1 = temp;

        Ng2 = std::max((int)std::sqrt(Np2/NpNg2),2);
        Nt2 = (double)Np2/NpNt;

        runtime1 = runtime2;
        runtime2 = Nt2*getRuntimePerStep((double)0.1/Nt2, Np2,Ng2);

        if(std::min(Np1,Np2) > 10000 && std::abs(Np1-Np2) < (double)std::max(Np1,Np2)/1000){
            break;
        }
        if(Np1 == Np2){
            break;
        }
    }

    return Np2 + (int)((totalRuntime-runtime2)*((double)(Np3-Np2))/(runtime3-runtime2));
}

int main(int argc, char* argv[]){
    int status = MPI_Init(&argc, &argv);
    if(status != MPI_SUCCESS){
        throw std::runtime_error("MPI initialisation error.");
    }

    std::vector<double> Powers = {-3,-2,-1,0,0.3,0.6,1,1.5,2,2.5,3,3.5,4};
    int id = 0;

    double totalRuntime = 0.1;
    double NpNg2 = 100;

    for(auto power: Powers){
        double NpNt = std::pow(10,power);
        int Np = searchNp(NpNt,NpNg2,totalRuntime);
        int Nt = (double)Np/NpNt;
        //int Ng = (int)std::sqrt(Np/NpNg2);
        int Ng = 30;

        if(Ng < 2){
            throw std::runtime_error("Not enough Ng");
        }

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)0.1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio2/semiimplicit1/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "Np.txt";
        std::ofstream npfile(ss2.str());
        npfile << Np;
        npfile.close();

        id += 1;
    }
    id = 0;

    Powers = {-2,-1,0,0.3,0.6,1,1.5,2,2.5,3,3.5,4};
    totalRuntime=0.01;
    for(auto power: Powers){
        double NpNt = std::pow(10,power);
        int Np = std::max(searchNp(NpNt,NpNg2,totalRuntime),2);
        int Nt = (double)Np/NpNt;
        //int Ng = (int)std::sqrt(Np/NpNg2);
        int Ng = 100;

        if(Ng < 2){
            throw std::runtime_error("Not enough Ng");
        }

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)0.1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio2/semiimplicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "Np.txt";
        std::ofstream npfile(ss2.str());
        npfile << Np;
        npfile.close();

        id += 1;
    }
    id = 0;

    Powers = {-3,-2,-1,0,0.3,0.6,1,1.5,2,2.5,3,3.5,4};

    totalRuntime = 0.1;
    NpNg2 = 100;

    for(auto power: Powers){
        double NpNt = std::pow(10,power);
        int Np = searchNp(NpNt,NpNg2,totalRuntime);
        int Nt = (double)Np/NpNt;
        //int Ng = (int)std::sqrt(Np/NpNg2);
        int Ng = 30;

        if(Ng < 2){
            throw std::runtime_error("Not enough Ng");
        }

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)0.1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio2/explicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "Np.txt";
        std::ofstream npfile(ss2.str());
        npfile << Np;
        npfile.close();

        id += 1;
    }
    id = 0;

    Powers = {-2,-1,0,0.3,0.6,1,1.5,2,2.5,3,3.5,4};
    totalRuntime=0.01;
    for(auto power: Powers){
        double NpNt = std::pow(10,power);
        int Np = std::max(searchNp(NpNt,NpNg2,totalRuntime),2);
        int Nt = (double)Np/NpNt;
        //int Ng = (int)std::sqrt(Np/NpNg2);
        int Ng = 100;

        if(Ng < 2){
            throw std::runtime_error("Not enough Ng");
        }

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)0.1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 0.1;
        std::stringstream ss;
        ss << "./outputs/optratio2/explicit2/" << id << "/";
        config.saveConfig.outputFilesDirectory = ss.str();
        Simulation<double,2,3> sim(config);
        sim.initialise();
        sim.run();

        std::stringstream ss2;
        ss2 << ss.str();
        ss2 << "Np.txt";
        std::ofstream npfile(ss2.str());
        npfile << Np;
        npfile.close();

        id += 1;
    }
    id = 0;


    MPI_Finalize();
}