#ifndef PIC_SEMI_IMPLICIT_PRESET_SIMULATIONS_H
#define PIC_SEMI_IMPLICIT_PRESET_SIMULATIONS_H

#include <vector>

#include "../models/Config.h"
#include "../models/Grid.h"
#include "math_helper.h"
#include "preset_distributions.h"
#include "preset_species.h"
#include "preset_fields.h"
#include "preset_save_configs.h"


namespace preset_configs {
    template <typename T>
    Config<T,1,1> landau1D1V();

    template <typename T>
    Config<T,2,3> constVelocityX();

    template <typename T>
    Config<T,2,3> landau2D3VX(unsigned int Nx, unsigned int Ny);

    //Same as previous case but initialised from file
    template <typename T>
    Config<T,2,3> landauFile();

    template <typename T>
    Config<T,2,3> landau2D3VXWave(unsigned int Nx, unsigned int Ny);

    template <typename T>
    Config<T,2,3> landau2D3VXWaveStatic(unsigned int Nx, unsigned int Ny);

    //Top hat species around 0,0 in a -1 to 1 domain which gets a constant electric acceleration.
    template <typename T>
    Config<T,2,3> constAccelerationX();

    //Constant magnetic field out of plane for constant gyration with radius 0.5
    template <typename T>
    Config<T,2,3> magneticGyration();

    //Same as previous but magnetic field is in X direction so rotation appears as oscillation in y.
    template <typename T>
    Config<T,2,3> magneticGyrationX();

    //Constant electric field towards the centre and proportional to the distance from the origin.
    template <typename T>
    Config<T,2,3> constPotentialWell();

    //Same as previous case but initialised from file
    template <typename T>
    Config<T,2,3> constPotentialWellFile();

    //Electron beam with no background ions
    template <typename T>
    Config<T,2,3> electronBeam(unsigned int Nx, unsigned int Ny);
}

#include "preset_configs.cpp"

#endif //PIC_SEMI_IMPLICIT_PRESET_SIMULATIONS_H
