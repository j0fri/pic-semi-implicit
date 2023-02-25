#include "Species.h"

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::Species(const Config<T,Nd,Nv>::SpeciesConfig& speciesConfig) : Np(speciesConfig.Np),
                                                                                   m(speciesConfig.m), q(speciesConfig.m), initialXDist(speciesConfig.xDist),
                                                                                   initialVDist(speciesConfig.vDist),
                                                                                   initialisePositionFromFile(speciesConfig.initialisePositionFromFile),
                                                                                   initialiseVelocityFromFile(speciesConfig.initialiseVelocityFromFile),
                                                                                   initialPositionFileName(speciesConfig.initialPositionFileName),
                                                                                   initialVelocityFileName(speciesConfig.initialVelocityFileName){}

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