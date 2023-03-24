//
// Created by jf1519 on 22/03/23.
//

#ifndef PIC_SEMI_IMPLICIT_FIELD2D3V_H
#define PIC_SEMI_IMPLICIT_FIELD2D3V_H

#include <vector>
#include <fstream>
#include "Config.h"
#include "Species.h"

template<typename T>
class Field2D3V {
private:
    T* field; /*Field is an array of length Nx*Ny*6, where Ex,Ey,Ez,Bx,By,Bz are stored successively for each grid point
    an increment of 6 is an increment in x value, and an increment of Nx is an increment in y value*/
    T* fieldT; //Field at half time step
    T* J; //Current is an array of length Nx*Ny*3, same structure as field with Jx,Jy,Jz

    T* Mgg; //Mass matrix for cell with itself
    T* Mgpg; //Mass matrix for cell with cell + dx
    T* Mggp; //Mass matrix for cell with cell + dy
    T* Mgpgp; //Mass matrix for cell with cell + dx + dy
    /*All the previous matrices are arrays of length Nx*Ny*9, same structure as field, with all 9 components of the
    matrix are stored in column-major format (increment of 1 is in x) for each grid value*/

    T* A; //System matrix: Nx*Ny*6*Nx*Ny*6
    T* c; //System vector: Nx*Ny*6, same order as field

public:
    Field2D3V() = delete;
    explicit Field2D3V(const Config<T,2,3>::FieldConfig& fieldConfig);
    Field2D3V(const Field2D3V& other);
    ~Field2D3V();

    void initialise(const std::vector<Species<T,2,3>*>& species);

    void saveElectricField(std::ofstream& outputFile) const ;
    void saveMagneticField(std::ofstream& outputFile) const;
    void saveEnergy(std::ofstream& outputFile) const;
    void saveVoltage(std::ofstream& outputFile) const;

    const T* getField() const;
    const T* getFieldT() const;

private:
    void accumulateJ(const std::vector<Species<T,2,3>*>& species);
    void accumulateM(const std::vector<Species<T,2,3>*>& species, T dt);
    void solveAndAdvance(T dt);
};


#endif //PIC_SEMI_IMPLICIT_FIELD2D3V_H
