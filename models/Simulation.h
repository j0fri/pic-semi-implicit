#ifndef PIC_SEMI_IMPLICIT_SIMULATION_H
#define PIC_SEMI_IMPLICIT_SIMULATION_H

#include "Config.h"
#include "Species.h"
#include "Species1D1V.h"
#include "Field.h"
#include "Field1D1V.h"


template <typename T, unsigned int Nd, unsigned int Nv>
class Simulation {
public:
    enum State{Uninitialised, Initialised, InitialisationError, Finalised, RuntimeError};
private:
    State state;
    std::vector<Species<T,Nd,Nv>*> species;
    Field<T,Nd,Nv>* field;
    Config<T,Nd,Nv>::TimeConfig timeConfig;
    Config<T,Nd,Nv>::SaveConfig saveConfig;
    bool verbose;

public:
    Simulation() = delete;
    explicit Simulation(const Config<T,Nd,Nv>& config);
    ~Simulation();
    void initialise();
    void clearOutputFiles();
    void run();
    void save();
    void checkValidState() const;
    State getState() const;
};


#endif //PIC_SEMI_IMPLICIT_SIMULATION_H
