#include <iostream>
#include "Simulation.h"


template <typename T, unsigned int Nd, unsigned int Nv>
Simulation<T,Nd,Nv>::Simulation(const Config<T,Nd,Nv>& config) {
    species = std::vector<Species<T,Nd,Nv>*>(config.speciesConfig.size());
    for(int i = 0; i < (int)species.size(); ++i){
        if(Nd == 1 && Nv == 1){
            species[i] = new Species1D1V<T>(config.speciesConfig[i], config.bcConfig);
//        }else if(Nd == 2 && Nv == 3){
//            species[i] = new Species2D3V<T>(config.species[i]);
        }else{
            state = State::InitialisationError;
            throw std::runtime_error("Nd and Nv combination not supported.");
        }
    }
    if(Nd == 1 && Nv == 1){
        field = new Field1D1V<T>(config.fieldConfig, config.bcConfig);
//    }else if(Nd == 2 && Nv == 3){
//        field = new Field2D3V<T>(config.fieldConfig);
    }else{
        state = State::InitialisationError;
        throw std::runtime_error("Nd and Nv combination not supported.");
    }
    state = State::Uninitialised;
    timeConfig = config.timeConfig;
    saveConfig = config.saveConfig;
    bcConfig = config.bcConfig;
    verbose = config.verbose;
}


template <typename T, unsigned int Nd, unsigned int Nv>
Simulation<T,Nd,Nv>::~Simulation(){
    delete field;
    for(Species<T,Nd,Nv>* sPtr: species){
        delete sPtr;
    }
}


template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T, Nd, Nv>::initialise() {
    try{
        for(int i = 0; i < (int)species.size(); ++i){
            species[i]->initialise();
        }
        field->initialise(species);
        this->clearOutputFiles();
        state = State::Initialised;
    }catch(const std::exception& exception){
        std::cout << "Exception thrown during initialisation: " << exception.what() << std::endl;
        state = State::InitialisationError;
    }
}


template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::clearOutputFiles(){
    if(saveConfig.savePosition){
        std::ofstream speciesPositionFile(saveConfig.outputFilesDirectory + saveConfig.speciesPositionFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!speciesPositionFile.is_open()){
            throw std::runtime_error("Could not open species position output file.");
        }
        speciesPositionFile.close();
    }

    if(saveConfig.savePositionDistribution){
        std::ofstream speciesPositionDistributionFile(saveConfig.outputFilesDirectory + saveConfig.speciesPositionDistributionFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!speciesPositionDistributionFile.is_open()){
            throw std::runtime_error("Could not open species position distribution output file.");
        }
        speciesPositionDistributionFile.close();
    }

    if(saveConfig.saveVelocity){
        std::ofstream speciesVelocityFile(saveConfig.outputFilesDirectory + saveConfig.speciesVelocityFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!speciesVelocityFile.is_open()){
            throw std::runtime_error("Could not open species velocity output file.");
        }
        speciesVelocityFile.close();
    }

    if(saveConfig.saveVelocityDistribution){
        std::ofstream speciesVelocityDistributionFile(saveConfig.outputFilesDirectory + saveConfig.speciesVelocityDistributionFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!speciesVelocityDistributionFile.is_open()){
            throw std::runtime_error("Could not open species velocity distribution output file.");
        }
        speciesVelocityDistributionFile.close();
    }
    
    if(saveConfig.saveSpeciesEnergy){
        std::ofstream speciesEnergyFile(saveConfig.outputFilesDirectory + saveConfig.speciesEnergyFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!speciesEnergyFile.is_open()){
            throw std::runtime_error("Could not open species energy file.");
        }
        speciesEnergyFile.close();
    }

    if(saveConfig.saveElectricField){
        std::ofstream electricFieldFile(saveConfig.outputFilesDirectory + saveConfig.electricFieldFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!electricFieldFile.is_open()){
            throw std::runtime_error("Could not open electric field file.");
        }
        electricFieldFile.close();
    }

    if(saveConfig.saveMagneticField){
        std::ofstream magneticFieldFile(saveConfig.outputFilesDirectory + saveConfig.magneticFieldFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!magneticFieldFile.is_open()){
            throw std::runtime_error("Could not open magnetic field file.");
        }
        magneticFieldFile.close();
    }
    
    if(saveConfig.saveFieldEnergy){
        std::ofstream fieldEnergyFile(saveConfig.outputFilesDirectory + saveConfig.fieldEnergyFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!fieldEnergyFile.is_open()){
            throw std::runtime_error("Could not open field energy file.");
        }
        fieldEnergyFile.close();
    }

    if(saveConfig.saveVoltage){
        std::ofstream voltageFile(saveConfig.outputFilesDirectory + saveConfig.voltageFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
        if(!voltageFile.is_open()){
            throw std::runtime_error("Could not open voltage file.");
        }
        voltageFile.close();
    }
}


//TODO: check no more calls needed
template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::run() {
    this->checkValidState();
    try{
        T t = 0;
        T nextSave = 0;
        while(t <= timeConfig.total){
            if(verbose){
                std::cout << "t: " << t << std::endl;
            }
            if (t >= nextSave){
                this->save();
                nextSave += saveConfig.saveInterval;
            }
            for(Species<T,Nd,Nv>* sPtr: species){
                sPtr->advancePositions(timeConfig.step, field);
            }
            field->advanceField(species, timeConfig.step);
            for(Species<T,Nd,Nv>* sPtr: species){
                sPtr->advanceVelocities(timeConfig.step, field);
            }
            t += timeConfig.step;
        }
    }catch(const std::exception& exception){
        std::cout << "Exception thrown during simulation: " << exception.what() << std::endl;
        state = State::RuntimeError;
    }
}


//TODO: add capability to save specific species only
template<typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::save(){
    if(saveConfig.savePosition){
        std::ofstream speciesPositionFile(saveConfig.outputFilesDirectory + saveConfig.speciesPositionFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!speciesPositionFile.is_open()){
            throw std::runtime_error("Could not open species position output file.");
        }
        for(unsigned int i = 0; i < (unsigned int)species.size(); ++i){
            species[i]->savePosition(speciesPositionFile);
        }
        speciesPositionFile.close();
    }

    if(saveConfig.savePositionDistribution){
        std::ofstream speciesPositionDistributionFile(saveConfig.outputFilesDirectory + saveConfig.speciesPositionDistributionFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!speciesPositionDistributionFile.is_open()){
            throw std::runtime_error("Could not open species position distribution output file.");
        }
        for(unsigned int i = 0; i < (unsigned int)species.size(); ++i){
            species[i]->savePositionDistribution(speciesPositionDistributionFile, field);
        }
        speciesPositionDistributionFile.close();
    }

    if(saveConfig.saveVelocity){
        std::ofstream speciesVelocityFile(saveConfig.outputFilesDirectory + saveConfig.speciesVelocityFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!speciesVelocityFile.is_open()){
            throw std::runtime_error("Could not open species velocity output file.");
        }
        for(unsigned int i = 0; i < (unsigned int)species.size(); ++i){
            species[i]->saveVelocity(speciesVelocityFile);
        }
        speciesVelocityFile.close();
    }

    if(saveConfig.saveVelocityDistribution){
        std::ofstream speciesVelocityDistributionFile(saveConfig.outputFilesDirectory + saveConfig.speciesVelocityDistributionFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!speciesVelocityDistributionFile.is_open()){
            throw std::runtime_error("Could not open species velocity distribution output file.");
        }
        for(unsigned int i = 0; i < (unsigned int)species.size(); ++i){
            species[i]->saveVelocityDistribution(speciesVelocityDistributionFile);
        }
        speciesVelocityDistributionFile.close();
    }

    if(saveConfig.saveSpeciesEnergy){
        std::ofstream speciesEnergyFile(saveConfig.outputFilesDirectory + saveConfig.speciesEnergyFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!speciesEnergyFile.is_open()){
            throw std::runtime_error("Could not open species energy file.");
        }
        for(unsigned int i = 0; i < (unsigned int)species.size(); ++i){
            species[i]->saveEnergy(speciesEnergyFile);
        }
        speciesEnergyFile.close();
    }

    if(saveConfig.saveElectricField){
        std::ofstream electricFieldFile(saveConfig.outputFilesDirectory + saveConfig.electricFieldFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!electricFieldFile.is_open()){
            throw std::runtime_error("Could not open electric field file.");
        }
        field->saveElectricField(electricFieldFile);
        electricFieldFile.close();
    }

    if(saveConfig.saveMagneticField){
        std::ofstream magneticFieldFile(saveConfig.outputFilesDirectory + saveConfig.magneticFieldFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!magneticFieldFile.is_open()){
            throw std::runtime_error("Could not open magnetic field file.");
        }
        field->saveMagneticField(magneticFieldFile);
        magneticFieldFile.close();
    }

    if(saveConfig.saveFieldEnergy){
        std::ofstream fieldEnergyFile(saveConfig.outputFilesDirectory + saveConfig.fieldEnergyFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!fieldEnergyFile.is_open()){
            throw std::runtime_error("Could not open field energy file.");
        }
        field->saveEnergy(fieldEnergyFile);
        fieldEnergyFile.close();
    }

    if(saveConfig.saveVoltage){
        std::ofstream voltageFile(saveConfig.outputFilesDirectory + saveConfig.voltageFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!voltageFile.is_open()){
            throw std::runtime_error("Could not open voltage file.");
        }
        field->saveVoltage(voltageFile);
        voltageFile.close();
    }
}


template<typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::checkValidState() const{
    if (state == State::Uninitialised){
        throw std::runtime_error("Simulation must be initialised before running.");
    }
    if (state == State::InitialisationError){
        throw std::runtime_error("Simulation cannot run as there was an initialisation error.");
    }
    if (state == State::Finalised){
        throw std::runtime_error("Simulation already ran.");
    }
    if (state == State::RuntimeError){
        throw std::runtime_error("Simulation already ran and exited with error.");
    }
}


template<typename T, unsigned int Nd, unsigned int Nv>
Simulation<T,Nd,Nv>::State Simulation<T,Nd,Nv>::getState() const {
    return state;
}


template class Simulation<double,1,1>;
template class Simulation<float,1,1>;