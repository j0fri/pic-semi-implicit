#ifndef PIC_SEMI_IMPLICIT_PRESET_SAVE_CONFIGS_H
#define PIC_SEMI_IMPLICIT_PRESET_SAVE_CONFIGS_H

#include "../models/Config.h"

namespace preset_save_configs{
    template<typename T, unsigned int Nd, unsigned int Nv>
    typename Config<T,Nd,Nv>::SaveConfig Energies(T dt);
}

#include "preset_save_configs.cpp"

#endif //PIC_SEMI_IMPLICIT_PRESET_SAVE_CONFIGS_H
