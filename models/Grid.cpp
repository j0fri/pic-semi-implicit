#include "Grid.h"

template <typename T, unsigned int Nd>
std::array<T,Nd> Grid<T,Nd>::getSpacings() const{
    std::array<T,Nd> out{};
    for(unsigned int i = 0; i < Nd; ++i){
        const auto& dim = dimensions[i];
        out[i] = (dim.max-dim.min)/dim.Nc;
    }
    return out;
}

template struct Grid<double,1>;
template struct Grid<double,2>;
template struct Grid<double,3>;
template struct Grid<float,1>;
template struct Grid<float,2>;
template struct Grid<float,3>;