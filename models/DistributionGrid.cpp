#include "DistributionGrid.h"

template<typename T, unsigned int Nd>
DistributionGrid<T, Nd>::DistributionGrid(const std::function<T(std::array<T, Nd>)> &f, const Grid<T, Nd> &grid):
        Distribution<T, Nd>(f), grid{grid} {}

template<typename T, unsigned int Nd>
DistributionGrid<T, Nd>::DistributionGrid(const Distribution<T, Nd> &dist, const Grid<T, Nd> &grid):
        Distribution<T, Nd>(dist.f), grid{grid} {}

template<typename T, unsigned int Nd>
std::vector<std::array<T, Nd>> DistributionGrid<T, Nd>::generate(int Np) const {
    return ((Distribution<T,Nd>*)this)->generate(Np, grid);
}

template class DistributionGrid<double,1>;
template class DistributionGrid<double,2>;
template class DistributionGrid<double,3>;
template class DistributionGrid<float,1>;
template class DistributionGrid<float,2>;
template class DistributionGrid<float,3>;