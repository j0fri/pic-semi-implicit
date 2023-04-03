#ifndef PIC_SEMI_IMPLICIT_OUTPUT_HELPER_H
#define PIC_SEMI_IMPLICIT_OUTPUT_HELPER_H

#include <iostream>

namespace output_helper {
    //Prints an m by n matrix stored in column-major format where ldx and incx are the leading dimension and increment
    //in the column-wise direction
    template <typename T>
    void outputColMajorMatrix(const T* a, unsigned int m, unsigned int n, unsigned int ldx, unsigned int incx,
                              std::ostream& ostream);

    //Prints an m by n matrix stored in row-major format where ldx and incx are the leading dimension and increment
    //in the rows-wise direction
    template <typename T>
    void outputRowMajorMatrix(const T* a, unsigned int m, unsigned int n, unsigned int ldx, unsigned int incx,
                              std::ostream& ostream);
}

#include "output_helper.cpp"

#endif //PIC_SEMI_IMPLICIT_OUTPUT_HELPER_H
