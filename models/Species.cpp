#include "Species.h"
#include <iostream>

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::Species(const typename Config<T,Nd,Nv>::SpeciesConfig& speciesConfig,
                            const typename Config<T,Nd,Nv>::BCConfig& bcConfig)
    : Np(Species::calculateLocalNp(speciesConfig.Np)), m(speciesConfig.m/speciesConfig.Np), q(speciesConfig.q/speciesConfig.Np),
    initialXDist(speciesConfig.xDist), initialXGrid(speciesConfig.initialXGrid), initialVDist(speciesConfig.vDist),
    initialVGrid(speciesConfig.initialVGrid), initialisePositionFromFile(speciesConfig.initialisePositionFromFile),
    initialPositionFileName(speciesConfig.initialPositionFileName), initialiseVelocityFromFile(speciesConfig.initialiseVelocityFromFile),
    initialVelocityFileName(speciesConfig.initialVelocityFileName), bcConfig(bcConfig),
    bcPositionGenerator(speciesConfig.bcPositionGenerator),
    bcNormalVelocityGenerator(speciesConfig.bcNormalVelocityGenerator), totalNp{speciesConfig.Np} {}

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::Species(const Species<T, Nd, Nv> &other) = default;

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::~Species() = default;

template<typename T, unsigned int Nd, unsigned int Nv>
void Species<T, Nd, Nv>::initialise(){
    if(initialisePositionFromFile){
        std::ifstream positionFile(initialPositionFileName, std::ios::in);
        if(!positionFile.is_open()){
            throw std::runtime_error("Could not open input position file.");
        }
        initialisePositions(positionFile);
        positionFile.close();
    }else{
        initialisePositions();
    }
    if(initialiseVelocityFromFile){
        std::ifstream velocityFile(initialVelocityFileName, std::ios::in);
        if(!velocityFile.is_open()){
            throw std::runtime_error("Could not open input velocity file.");
        }
        initialiseVelocities(velocityFile);
        velocityFile.close();
    }else{
        initialiseVelocities();
    }
}

template<typename T, unsigned int Nd, unsigned int Nv>
unsigned int Species<T, Nd, Nv>::calculateLocalNp(unsigned int Np) {
    int processId, numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    if(numProcesses>Np){
        throw std::invalid_argument("Number of particles must be larger than the number of processes.");
    }

    unsigned int baseNp = Np/numProcesses;
    unsigned int missing = Np-numProcesses*baseNp;
    if(processId < missing){
        return baseNp + 1;
    }else{
        return baseNp;
    }
}


template class Species<float,1,1>;
//template class Species<float,1,3>;
template class Species<float,2,3>;
template class Species<double,1,1>;
//template class Species<double,1,3>;
template class Species<double,2,3>;