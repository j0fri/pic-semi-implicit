#include "Simulation.h"


template <typename T, unsigned int Nd, unsigned int Nv>
Simulation<T,Nd,Nv>::Simulation(const Config<T,Nd,Nv>& config) {
    species = std::vector<Species<T,Nd,Nv>*>(config.speciesConfig.size());
    for(int i = 0; i < (int)species.size(); ++i){
        if(Nd == 1 && Nv == 1){
            species[i] = new Species1D1V(config.species[i]);
//        }else if(Nd == 2 && Nv == 3){
//            species[i] = new Species2D3V(config.species[i]);
        }else{
            state = State::InitialisationError;
            throw std::runtime_error("Nd and Nv combination not supported.");
        }
    }
    if(Nd == 1 && Nv == 1){
        field = new Field1D1V(config.fieldConfig);
//    }else if(Nd == 2 && Nv == 3){
//        field = new Field2D3V(config.fieldConfig);
    }else{
        state = State::InitialisationError;
        throw std::runtime_error("Nd and Nv combination not supported.");
    }
    state = State::Uninitialised;
    timeConfig = config.timeConfig;
    saveConfig = config.saveConfig;
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
        state = State::InitialisationError;
    }
}


template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::clearOutputFiles(){
    std::ofstream speciesOutputFile(saveConfig.speciesOutputFileName, std::ios::trunc);
    if(!speciesOutputFile.is_open()){
        throw std::runtime_error("Could not open species output file.");
    }
    speciesOutputFile.close();
    std::ofstream fieldOutputFile(saveConfig.fieldOutputFileName, std::ios::trunc);
    if(!speciesOutputFile.is_open()){
        throw std::runtime_error("Could not open field output file.");
    }
    fieldOutputFile.close();
}


//TODO: check no more calls needed
template <typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::run() {
    this->checkValidState();
    try{
        T t = 0;
        T nextSave = 0;
        while(t <= timeConfig.total){
            if (t >= nextSave){
                this->save();
                nextSave += saveConfig.saveInterval;
            }
            for(Species<T,Nd,Nv>* sPtr: species){
                sPtr->advancePositions(timeConfig.step, field);
            }
            field->advanceField(species);
            for(Species<T,Nd,Nv>* sPtr: species){
                sPtr->advanceVelocities(timeConfig.step, field);
            }
            t += timeConfig.step;
        }
    }catch(const std::exception& exception){
        state = State::RuntimeError;
    }
}


//TODO: add capability to save specific species only
template<typename T, unsigned int Nd, unsigned int Nv>
void Simulation<T,Nd,Nv>::save(){
    std::ofstream speciesOutputFile("outputs/species.txt",std::ios::app);
    if(!speciesOutputFile.is_open()){
        throw std::runtime_error("Could not open species output file.");
    }
    for(int i = 0; i < (int)species.size(); ++i){
        species[i]->save(speciesOutputFile,saveConfig);
    }
    speciesOutputFile.close();

    std::ofstream fieldOutputFile(saveConfig.fieldOutputFileName, std::ios::app);
    if(!speciesOutputFile.is_open()){
        throw std::runtime_error("Could not open field output file.");
    }
    field->save(fieldOutputFile,saveConfig);
    fieldOutputFile.close();
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