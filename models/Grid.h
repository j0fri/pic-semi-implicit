#ifndef PIC_SEMI_IMPLICIT_GRID_H
#define PIC_SEMI_IMPLICIT_GRID_H

#include<array>

template <typename T, unsigned int Nd>
struct Grid {
    struct Dim{
        T min;
        T max;
        unsigned int Nc;
    };
    std::array<Dim,Nd> dimensions;
    std::array<T,Nd> getSpacings() const;
};


#endif //PIC_SEMI_IMPLICIT_GRID_H
