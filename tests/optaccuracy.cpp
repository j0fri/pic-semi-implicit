#include <iostream>

#include "../models/Config.h"
#include "../models/Simulation.h"
#include "../helpers/preset_configs.h"
#include "../helpers/output_helper.h"

double getRuntimePerStep(double dt, int Np, int Ng, bool exp = false){
    std::cout << "Np=" << Np << std::endl;
    auto config = preset_configs::langmuir(Np,Ng/2,2,dt);
    config.saveConfig.outputFilesDirectory = "./outputs/optaccuracy/dummy/";
    config.saveConfig.saveAllTimes = true;
    config.timeConfig.total = dt*100;
    config.useExplicitScheme = exp;

    Simulation<double,2,3> sim(config);
    sim.initialise();
    sim.run();

    std::ifstream dummyRuntimes("./outputs/optaccuracy/dummy/runtime.txt");
    if(!dummyRuntimes.is_open()){
        throw std::runtime_error("Dummy runtime file not open.");
    }

    double initTime = 0;
    double runTime = 0;
    dummyRuntimes >> initTime;
    dummyRuntimes >> runTime;
    dummyRuntimes.close();

    return runTime/100;
}

int searchNp(double NpNt, double NpNg2, double totalRuntime){
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

    int testNp = 100000;
//    std::vector<double> Nts = {1,2,5,10,30,60,100,250,500,1000};
//    std::vector<double> Nts = {5,7,10,15,22,30,45,60,80,100,250,500,750,1000,1250,1500,2000,2500};
    std::vector<double> Runtimes = {0.000001,0.000003,0.00001,0.00003,0.0001,0.003,0.000001,0.000003,0.00001,0.00003,0.0001,0.003,0.000001,0.000003,0.00001,0.00003,0.0001,0.003,0.001,0.003,0.01,0.03,0.001,0.003,0.01,0.03,0.001,0.003,0.01,0.03,0.1,0.2,0.3,0.5,1};
    int id = 0;

    double runtimePerStep = getRuntimePerStep(0.01,testNp,30);
    double C = testNp/runtimePerStep;

    double NpNg2 = 100;
    double NpNt = 0.1;

    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)std::sqrt((double)Np/NpNg2),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/semiimplicit1/" << id << "/";
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

    NpNt = 1;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)std::sqrt((double)Np/NpNg2),30);
        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/semiimplicit2/" << id << "/";
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

    NpNt = 10;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)std::sqrt((double)Np/NpNg2),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/semiimplicit3/" << id << "/";
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

    NpNt = 100;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)std::sqrt((double)Np/NpNg2),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/semiimplicit4/" << id << "/";
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

    NpNt = 1000;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)std::sqrt((double)Np/NpNg2),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/semiimplicit5/" << id << "/";
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



    runtimePerStep = getRuntimePerStep(0.01,testNp,30, true);
    C = testNp/runtimePerStep;
    double NpNg = 10000;
    NpNt = 0.1;

    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)((double)Np/NpNg),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/explicit1/" << id << "/";
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

    NpNt = 1;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)((double)Np/NpNg),30);
        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/explicit2/" << id << "/";
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

    NpNt = 10;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)((double)Np/NpNg),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/explicit3/" << id << "/";
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

    NpNt = 100;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)((double)Np/NpNg),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/explicit4/" << id << "/";
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

    NpNt = 1000;
    for(auto runtime: Runtimes){
        int Np = (int)std::sqrt(C*runtime*NpNt);
        int Nt = std::max((int) ((double)Np/NpNt),1);
        if(Np < 5){
//            throw std::runtime_error("Not enough/too many Np");
            continue;
        }
        if(Np > 1000000){
            continue;
        }
        int Ng = std::max((int)((double)Np/NpNg),30);

        auto config = preset_configs::langmuir(Np,Ng/2,2,(double)1/Nt);
        config.saveConfig.saveAllTimes = true;
        config.timeConfig.total = 1;
        config.useExplicitScheme = true;
        std::stringstream ss;
        ss << "./outputs/optaccuracy/explicit5/" << id << "/";
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