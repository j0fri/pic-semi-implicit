#include "Field.h"
#include "Config.h"

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::Field(const Config<T,Nd,Nv>::FieldConfig& fieldConfig) : grid(fieldConfig.grid), c(fieldConfig.c),
                e0(fieldConfig.e0), forcedE(fieldConfig.forcedE), forcedB(fieldConfig.forcedB),
                onlyForcedE(fieldConfig.onlyForcedE), onlyForcedB(fieldConfig.onlyForcedB),
                initialiseFromSpecies(fieldConfig.initialiseFromSpecies){}

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::Field(const Field<T,Nd,Nv>& other) = default;

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::~Field() = default;

template<typename T, unsigned int Nd, unsigned int Nv>
void Field<T, Nd, Nv>::advanceField(const std::vector<Species<T,Nd,Nv>*>& species, T dt) {
    this->accumulateJ(species);
    this->accumulateM(species, dt);
    this->solveAndAdvance(dt);
}

template class Field<float,1,1>;
//template class Field<float,1,3>;
//template class Field<float,2,3>;
template class Field<double,1,1>;
//template class Field<double,1,3>;
//template class Field<double,2,3>;