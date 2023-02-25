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