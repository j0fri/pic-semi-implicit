#ifndef PIC_SEMI_IMPLICIT_SIMULATION_H
#define PIC_SEMI_IMPLICIT_SIMULATION_H

#include "Config.h"
#include "Species.h"
#include "Species1D1V.h"
#include "Species2D3V.h"
#include "Field.h"
#include "Field1D1V.h"
#include "Field2D3V.h"
#include "Field2D3VExplicit.h"


template <typename T, unsigned int Nd, unsigned int Nv>
class Simulation {
public:
    enum State{Uninitialised, Initialised, InitialisationError, Finalised, RuntimeError};
private:
    int processId;
    int numProcesses;
    State state;
    std::vector<Species<T,Nd,Nv>*> species;
    Field<T,Nd,Nv>* field;
    typename Config<T,Nd,Nv>::TimeConfig timeConfig;
    typename Config<T,Nd,Nv>::SaveConfig saveConfig;
    typename Config<T,Nd,Nv>::BCConfig bcConfig;
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

private:
    void outputConfig(const Config<T,Nd,Nv>& config) const;
};


#endif //PIC_SEMI_IMPLICIT_SIMULATION_H
