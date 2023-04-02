#include "output_helper.h"

template <typename T>
void output_helper::outputColMajorMatrix(const T* a, int m, int n, int ldx, int incx, std::ostream& ostream){
    for(unsigned int row = 0; row < m; ++row){
        for(unsigned int col = 0; col < n; ++col){
            ostream << a[incx*row+ldx*col];
            if(col != n-1){
                ostream << " ";
            }
        }
        ostream << std::endl;
    }
}

template <typename T>
void output_helper::outputRowMajorMatrix(const T* a, int m, int n, int ldx, int incx, std::ostream& ostream){
    for(unsigned int row = 0; row < m; ++row){
        for(unsigned int col = 0; col < n; ++col){
            ostream << a[incx*col+ldx*row];
            if(col != n-1){
                ostream << " ";
            }
        }
        ostream << std::endl;
    }
}
