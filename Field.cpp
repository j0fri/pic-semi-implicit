#include <Field.h>

template <typename T>
Field<T>::Field(int Nx, T Lx, T c, T e0): Nd(0), Nx(Nx), Lx(Lx), dx(Lx/Nx), c(c), e0(e0) {};

template <typename T>
Field<T>::Field(const Field& other): Nd(other.Nd), Nx(other.Nx), Lx(other.Lx), dx(other.dx), c(other.c), e0(other.e0) {};

template <typename T>
Field<T>::~Field() {};