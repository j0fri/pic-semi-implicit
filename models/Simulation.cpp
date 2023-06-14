#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <string>

#include <chrono>
#include <filesystem>

#include "Simulation.h"
#include "../helpers/string_helper.h"


template <typename T, unsigned int Nd, unsigned int Nv>
Simulation<T,Nd,Nv>::Simulation(const Config<T,Nd,Nv>& config): processId(0), numProcesses(0) {
    try{
        species = std::vector<Species<T,Nd,Nv>*>(config.speciesConfig.size());
        for(int i = 0; i < (int)species.size(); ++i){
            if constexpr (Nd == 1 && Nv == 1) {
                species[i] = new Species1D1V<T>(config.speciesConfig[i], config.bcConfig);
            }else if constexpr(Nd == 2 && Nv == 3){
                species[i] = new Species2D3V<T>(config.speciesConfig[i], config.bcConfig);
            }else{
                state = State::InitialisationError;
                throw std::invalid_argument("Nd and Nv combination not supported.");
            }
        }
    }catch(const std::bad_alloc& e){
        std::cout << "Memory allocation exception during species construction." << std::endl;
        throw;
    }

    try{
        if constexpr (Nd == 1 && Nv == 1) {
            if(!config.useExplicitScheme){
                field = new Field1D1V<T>(config.fieldConfig, config.bcConfig);
            }else{
                throw std::invalid_argument("Explicit scheme not implemented in 1D1V.");
            }
        } else if constexpr (Nd == 2 && Nv == 3){

                if(!config.useExplicitScheme){
                    field = new Field2D3V<T>(config.fieldConfig, config.bcConfig);
                }else{
                    field = new Field2D3VExplicit<T>(config.fieldConfig, config.bcConfig);
                }
        }else{
            state = State::InitialisationError;
            throw std::runtime_error("Nd and Nv combination not supported.");
        }
    }catch(const std::bad_alloc& e){
        std::cout << "Memory allocation exception during field construction." << std::endl;
        throw;
    }

    state = State::Uninitialised;
    timeConfig = config.timeConfig;
    saveConfig = config.saveConfig;
    bcConfig = config.bcConfig;
    verbose = config.verbose;

    //Create output folder if it doesn't exist
    std::filesystem::create_directories(saveConfig.outputFilesDirectory);
    if(config.outputConfig){
        this->outputConfig(config);
    }
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
    auto begin = std::chrono::steady_clock::now();
    try{
        for(int i = 0; i < (int)species.size(); ++i){
            species[i]->initialise();
            species[i]->advancePositions(0, field); //Important for weight calculations
        }
        field->initialise(species, timeConfig.step);
        this->clearOutputFiles();
        state = State::Initialised;
    }catch(const std::exception& exception){
        std::cout << "Exception thrown during initialisation: " << exception.what() << std::endl;
        state = State::InitialisationError;
    }
    auto end = std::chrono::steady_clock::now();
    initialisationTime = (T)std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/10e9;
}


template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::clearOutputFiles(){
    if(processId == 0){
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

        if(saveConfig.saveElectrostaticPotential){
            std::ofstream electrostaticPotentialFile(saveConfig.outputFilesDirectory + saveConfig.electrostaticPotentialFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
            if(!electrostaticPotentialFile.is_open()){
                throw std::runtime_error("Could not open electrostatic potential file.");
            }
            electrostaticPotentialFile.close();
        }

        if(saveConfig.saveCurrent){
            std::ofstream currentFile(saveConfig.outputFilesDirectory + saveConfig.currentFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
            if(!currentFile.is_open()){
                throw std::runtime_error("Could not open current file.");
            }
            currentFile.close();
        }

        if(saveConfig.saveRuntime){
            std::ofstream runtimeFile(saveConfig.outputFilesDirectory + saveConfig.runtimeFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
            if(!runtimeFile.is_open()){
                throw std::runtime_error("Could not open runtime file.");
            }
            runtimeFile.close();
        }
        
        if(saveConfig.saveSolverSteps){
            std::ofstream solverStepsFile(saveConfig.outputFilesDirectory + saveConfig.solverStepsFileName + saveConfig.outputFilesSubscript, std::ios::trunc);
            if(!solverStepsFile.is_open()){
                throw std::runtime_error("Could not open solverSteps file.");
            }
            solverStepsFile.close();
        }
    }
}


template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::run() {
    auto begin = std::chrono::steady_clock::now();
    auto tempCounter = begin;
    auto tempCounter2 = begin;
    systemSolveTime = 0;
    this->checkValidState();
    try{
        T t = 0;
        T nextSave = 0;
        while(t <= timeConfig.total){
            if(verbose && processId == 0){
                std::cout << "t: " << t << std::endl;
            }
            if (saveConfig.saveAllTimes || t >= nextSave){
                this->save();
                nextSave += saveConfig.saveInterval;
            }
            for(Species<T,Nd,Nv>* sPtr: species){
                sPtr->advancePositions(timeConfig.step, field);
            }
            field->accumulateParticles(species, timeConfig.step);

            tempCounter = std::chrono::steady_clock::now();
            field->advanceField(timeConfig.step);
            tempCounter2 = std::chrono::steady_clock::now();
            systemSolveTime += (T)std::chrono::duration_cast<std::chrono::nanoseconds>(tempCounter2-tempCounter).count()/10e9;

            for(Species<T,Nd,Nv>* sPtr: species){
                sPtr->advanceVelocities(timeConfig.step, field);
            }
            t += timeConfig.step;
        }
        state = State::Finalised;
    }catch(const std::exception& exception){
        std::cout << "Exception thrown during simulation: " << exception.what() << std::endl;
        state = State::RuntimeError;
    }
    auto end = std::chrono::steady_clock::now();
    simulationTime = (T)std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/10e9;
    if(state == State::Finalised){
        this->saveRuntime();
    }
}


//TODO: add capability to save specific species only
template<typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::save(){
    if(numProcesses > 1){
        if(saveConfig.savePosition){
            throw std::invalid_argument("savePosition not allowed with more than one process");
        }
        if(saveConfig.savePositionDistribution){
            throw std::invalid_argument("savePositionDistribution not allowed with more than one process");
        }
        if(saveConfig.saveVelocity){
            throw std::invalid_argument("saveVelocity not allowed with more than one process");
        }
        if(saveConfig.saveVelocityDistribution){
            throw std::invalid_argument("saveVelocityDistribution not allowed with more than one process");
        }
    }

    if(processId == 0){
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
        if(saveConfig.saveCurrent){
            std::ofstream currentFile(saveConfig.outputFilesDirectory + saveConfig.currentFileName + saveConfig.outputFilesSubscript, std::ios::app);
            if(!currentFile.is_open()){
                throw std::runtime_error("Could not open current file.");
            }
            field->saveCurrent(currentFile);
            currentFile.close();
        }
        if(saveConfig.saveSolverSteps){
            std::ofstream solverStepsFile(saveConfig.outputFilesDirectory + saveConfig.solverStepsFileName + saveConfig.outputFilesSubscript, std::ios::app);
            if(!solverStepsFile.is_open()){
                throw std::runtime_error("Could not open solverSteps file.");
            }
            int lastSolverSteps = field->getSolverSteps();
            solverStepsFile << lastSolverSteps << std::endl;
            solverStepsFile.close();
        }
    }
    //Save electrostatic potential runs on all processes as all the charges are needed.
    if(saveConfig.saveElectrostaticPotential){
        std::ofstream electrostaticPotentialFile(saveConfig.outputFilesDirectory + saveConfig.electrostaticPotentialFileName + saveConfig.outputFilesSubscript, std::ios::app);
        if(!electrostaticPotentialFile.is_open()){
            throw std::runtime_error("Could not open electrostatic potential file.");
        }
        field->saveElectrostaticPotential(electrostaticPotentialFile, species);
        electrostaticPotentialFile.close();
    }
    //Save species energy runs on all processes as all the particles are needed.
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
typename Simulation<T,Nd,Nv>::State Simulation<T,Nd,Nv>::getState() const {
    return state;
}

template<typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T, Nd, Nv>::outputConfig(const Config<T, Nd, Nv> &config) const {
    if constexpr (!(Nd==2&&Nv==3)){
        throw std::invalid_argument("Config output is only available for 2D3V simulations.");
    }
    std::ofstream configFile(saveConfig.outputFilesDirectory + "config" + saveConfig.outputFilesSubscript, std::ios::trunc);
    if(!configFile.is_open()){
        throw std::runtime_error("Could not open config file.");
    }

    std::vector<std::array<std::string,2>> saveQueue{};
    saveQueue.emplace_back(std::array<std::string,2>{"T", string_helper::toString(config.timeConfig.total)});
    saveQueue.emplace_back(std::array<std::string,2>{"dt", string_helper::toString(config.timeConfig.step)});
    saveQueue.emplace_back(std::array<std::string,2>{"saveInterval", string_helper::toString(config.saveConfig.saveInterval)});
    saveQueue.emplace_back(std::array<std::string,2>{"xmin", string_helper::toString(config.fieldConfig.grid.dimensions[0].min)});
    saveQueue.emplace_back(std::array<std::string,2>{"xmax", string_helper::toString(config.fieldConfig.grid.dimensions[0].max)});
    saveQueue.emplace_back(std::array<std::string,2>{"Nx", string_helper::toString(config.fieldConfig.grid.dimensions[0].Nc)});
    saveQueue.emplace_back(std::array<std::string,2>{"ymin", string_helper::toString(config.fieldConfig.grid.dimensions[1].min)});
    saveQueue.emplace_back(std::array<std::string,2>{"ymax", string_helper::toString(config.fieldConfig.grid.dimensions[1].max)});
    saveQueue.emplace_back(std::array<std::string,2>{"Ny", string_helper::toString(config.fieldConfig.grid.dimensions[1].Nc)});
    saveQueue.emplace_back(std::array<std::string,2>{"Ns", string_helper::toString((unsigned int)config.speciesConfig.size())});

    for(unsigned int row = 0; row < 2; ++row){
        for(unsigned int i = 0; i < (unsigned int)saveQueue.size(); ++i){
            configFile << saveQueue[i][row];
            if(i < (unsigned int)saveQueue.size() - 1){
                configFile << " ";
            }
        }
        if(row < 1){
            configFile << std::endl;
        }
    }
}

template<typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T, Nd, Nv>::saveRuntime() {
    if(processId == 0){
        if(saveConfig.saveRuntime){
            std::ofstream runtimeFile(saveConfig.outputFilesDirectory + saveConfig.runtimeFileName + saveConfig.outputFilesSubscript, std::ios::app);
            if(!runtimeFile.is_open()){
                throw std::runtime_error("Could not open runtime file.");
            }
            runtimeFile << std::setprecision(10) << initialisationTime << " " << simulationTime << " " << systemSolveTime << std::endl;
            runtimeFile.close();
        }
    }
}

template class Simulation<double,1,1>;
template class Simulation<float,1,1>;
template class Simulation<double,2,3>;
template class Simulation<float,2,3>;