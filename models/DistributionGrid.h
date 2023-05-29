#ifndef PIC_SEMI_IMPLICIT_DISTRIBUTIONGRID_H
#define PIC_SEMI_IMPLICIT_DISTRIBUTIONGRID_H

#include "Distribution.h"

template <typename T, unsigned int Nd>
class DistributionGrid: Distribution<T,Nd> {
public:
    Grid<T,Nd> grid;
    explicit DistributionGrid(const std::function<T(std::array<T,Nd>)>& f, const Grid<T,Nd>& grid);
    explicit DistributionGrid(const Distribution<T,Nd>& dist, const Grid<T,Nd>& grid);
    std::vector<std::array<T, Nd>> generate(int Np) const;
};


#endif //PIC_SEMI_IMPLICIT_DISTRIBUTIONGRID_H
