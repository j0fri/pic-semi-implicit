#ifndef PIC_SEMI_IMPLICIT_PRESET_DISTRIBUTIONS_H
#define PIC_SEMI_IMPLICIT_PRESET_DISTRIBUTIONS_H

#include "../models/Distribution.h"
#include "../models/DistributionGrid.h"

namespace preset_distributions{
    //Boltzmann distribution in Nd dimensions, inputs are mass, boltzmann constant and temperature
    template <typename T, unsigned int Nd>
    Distribution<T,Nd> Boltzmann(T m, T Kb, T T0);

    //Uniform distribution in Nd dimension, input is the constant
    template <typename T, unsigned int Nd>
    Distribution<T,Nd> Uniform(T a);

    //Step distribution, first value is index of coordinate (1 for x-step, 2 for y-step...)
    //second value is the coordinate of the step value, and third is bool which if true gives
    //a 1 to 0 step instead of 0 to 1.
    template <typename T, unsigned int Nd>
    Distribution<T,Nd> Step(unsigned int Id, T a, bool reverse);

    //Step distribution, first value is index of coordinate (1 for x-sin, 2 for y-sin...)
    //rest of values are parameters in a*sin(u*b+c) where u=x,y,.. depending on Id
    template <typename T, unsigned int Nd>
    Distribution<T,Nd> Sin(unsigned int Id, T a, T b, T c);

    template <typename T, unsigned int Nd>
    Distribution<T,Nd> TopHat(unsigned int Id, T x1, T x2);

    //Constant value along domain
    template <typename T, unsigned int Nd>
    DistributionGrid<T,Nd> Constant(T val);
}

#include "preset_distributions.cpp"

#endif //PIC_SEMI_IMPLICIT_PRESET_DISTRIBUTIONS_H
