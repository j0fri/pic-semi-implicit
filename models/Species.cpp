#include "Species.h"

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::Species(unsigned int Np, T m, T q) : Np(Np), m(m), q(q){}

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::Species(const Species<T, Nd, Nv> &other) = default;

template<typename T, unsigned int Nd, unsigned int Nv>
Species<T, Nd, Nv>::~Species() = default;