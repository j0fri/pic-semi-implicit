#include "Species.h"

template <typename T>
Species<T>::Species(int Np, T m, T q): Nd(0), Nv(0), Np(Np), m(m), q(q) {};

template <typename T>
Species<T>::Species(const Species<T>& other): Nd(other.Nd), Nv(other.Nv), Np(other.Np), m(other.m), q(other.q) {};

template <typename T>
Species<T>::~Species() {};