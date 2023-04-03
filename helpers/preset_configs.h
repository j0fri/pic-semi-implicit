#ifndef PIC_SEMI_IMPLICIT_PRESET_SIMULATIONS_H
#define PIC_SEMI_IMPLICIT_PRESET_SIMULATIONS_H

#include <vector>

#include "../models/Config.h"
#include "../models/Grid.h"
#include "math_helper.h"
#include "preset_distributions.h"
#include "preset_species.h"
#include "preset_fields.h"


namespace preset_configs {
    template <typename T>
    Config<T,1,1> landau1D1V();

    template <typename T>
    Config<T,2,3> landau2D3VX(unsigned int Nx, unsigned int Ny);

    //Top hat species around 0,0 in a -1 to 1 domain which gets a constant electric acceleration.
    template <typename T>
    Config<T,2,3> constAccelerationX();

    template <typename T>
    Config<T,2,3> magneticGyration();
}

#include "preset_configs.cpp"

#endif //PIC_SEMI_IMPLICIT_PRESET_SIMULATIONS_H
