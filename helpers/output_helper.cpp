#include "output_helper.h"

#ifndef PIC_SEMI_IMPLICIT_OUTPUT_HELPER_CPP
#define PIC_SEMI_IMPLICIT_OUTPUT_HELPER_CPP

template <typename T>
void output_helper::outputColMajorMatrix(const T* a, unsigned int m, unsigned int n, unsigned int ldx,
                                         unsigned int incx, std::ostream& ostream){
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
void output_helper::outputRowMajorMatrix(const T* a, unsigned int m, unsigned int n, unsigned int ldx,
                                         unsigned int incx, std::ostream& ostream){
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

template <typename T>
bool output_helper::testSameFileContent(std::ifstream &file1, std::ifstream &file2, T tolerance) {
    if(!file1.is_open() || !file2.is_open()){
        throw std::invalid_argument("Compared files not open.");
    }
    T d1;
    T d2;
    while(true){
        if(file1.eof() && file2.eof()){
            return true;
        }
        if(file1.eof() || file2.eof()){
            return false;
        }
        file1 >> d1;
        file2 >> d2;
        if(std::abs(d1 - d2) > tolerance){
            return false;
        }
    }
}

#endif //PIC_SEMI_IMPLICIT_OUTPUT_HELPER_CPP