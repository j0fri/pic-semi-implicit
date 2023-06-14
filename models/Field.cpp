#include <mpi.h>
#include "Field.h"
#include "Config.h"

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::Field(const typename Config<T,Nd,Nv>::FieldConfig& fieldConfig,
                      const typename Config<T,Nd,Nv>::BCConfig& bcConfig)
    : processId(0), numProcesses(0), grid(fieldConfig.grid), c(fieldConfig.c),
    e0(fieldConfig.e0), forcedE(fieldConfig.forcedE), forcedB(fieldConfig.forcedB),
    onlyForcedE(fieldConfig.onlyForcedE), onlyForcedB(fieldConfig.onlyForcedB),
    initialiseFromSpecies(fieldConfig.initialiseFromSpecies), bcConfig(bcConfig),
    solverTolerance(fieldConfig.solverTolerance.has_value() ?  fieldConfig.solverTolerance.value() : (T)1e-8){

    int rankStatus = MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    int sizeStatus = MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    if(rankStatus != MPI_SUCCESS || sizeStatus != MPI_SUCCESS){
        throw std::runtime_error("Could not obtain MPI rank during field construction.");
    }
}

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::Field(const Field<T,Nd,Nv>& other) = default;

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::~Field() = default;

template<typename T, unsigned int Nd, unsigned int Nv>
void Field<T, Nd, Nv>::accumulateParticles(const std::vector<Species<T, Nd, Nv> *> &species, T dt) {
    this->accumulateJ(species);
    this->accumulateM(species, dt);
}

template<typename T, unsigned int Nd, unsigned int Nv>
void Field<T, Nd, Nv>::advanceField(T dt) {
    this->joinProcesses();
    this->solveAndAdvance(dt);
    this->distributeProcesses();
}

template<typename T, unsigned int Nd, unsigned int Nv>
int Field<T, Nd, Nv>::getSolverSteps() const {
    throw std::invalid_argument("Save solver steps not allowed for selected field type.");
}

template class Field<float,1,1>;
//template class Field<float,1,3>;
template class Field<float,2,3>;
template class Field<double,1,1>;
//template class Field<double,1,3>;
template class Field<double,2,3>;