#include "Field.h"

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::Field(const Grid<T,Nd>& grid, T c, T e0) : grid(grid), c(c), e0(e0) {}

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::Field(const Field<T,Nd,Nv>& other) = default;

template <typename T, unsigned int Nd, unsigned int Nv>
Field<T,Nd,Nv>::~Field() = default;