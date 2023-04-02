#ifndef PIC_SEMI_IMPLICIT_PRESET_FIELDS_H
#define PIC_SEMI_IMPLICIT_PRESET_FIELDS_H

#include "../models/Config.h"

namespace preset_fields {
    //Default1D1V is initialised from species charge distribution and has no forced/added term
    template <typename T>
    typename Config<T,1,1>::FieldConfig Default1D1V(T Lx, unsigned int Nx, T c, T e0);

    //Default2D3V is initialised from species charge distribution and has no forced/added term
    template <typename T>
    typename Config<T,2,3>::FieldConfig Default2D3V(T Lx, T Ly, unsigned int Nx, unsigned int Ny, T c, T e0);
}

#include "preset_fields.cpp"

#endif //PIC_SEMI_IMPLICIT_PRESET_FIELDS_H
