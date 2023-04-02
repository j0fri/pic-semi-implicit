#ifndef PIC_SEMI_IMPLICIT_PRESET_SPECIES_H
#define PIC_SEMI_IMPLICIT_PRESET_SPECIES_H

#include "../models/Config.h"

namespace preset_species {
    template <typename T>
    typename Config<T,1,1>::SpeciesConfig Uniform1D1V(unsigned int Np, T m, T q, T Lx, unsigned int Nx, T Kb, T T0);

    template <typename T>
    typename Config<T,2,3>::SpeciesConfig Uniform2D3V(unsigned int Np, T m, T q, T Lx, T Ly, unsigned int Nx, unsigned int Ny, T Kb, T T0);
};

#include "preset_species.cpp"

#endif //PIC_SEMI_IMPLICIT_PRESET_SPECIES_H
