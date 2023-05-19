#ifndef PIC_SEMI_IMPLICIT_FIELD2D3V_H
#define PIC_SEMI_IMPLICIT_FIELD2D3V_H

#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeLinearSolvers>

#include "Field.h"
#include "Config.h"
#include "Species.h"

template<typename T>
class Field2D3V: public Field<T,2,3> {
    typedef Eigen::SparseMatrix<T> SpMat;
public:
    const unsigned int Ng;
    const unsigned int Nx;
    const unsigned int Ny;
protected:
    T* field; /*Field is an array of length Nx*Ny*6, where Ex,Ey,Ez,Bx,By,Bz are stored successively for each grid point
    an increment of 6 is an increment in x value, and an increment of 6*Nx is an increment in y value. Electric field
    is stored at grid values while the magnetic field is stored at grid values shifted up and right by dx/2 and dy/2. */
    T* fieldT; //Field at half time step
    T* J; //Current is an array of length Nx*Ny*3, same structure as field with Jx,Jy,Jz

    //TODO: mass matrices and current are assumed to be centred at the electric field locations, not magnetic
    T* Mg; //Mass matrix for cell with itself
    T* Mgdx; //Mass matrix for cell with cell + dx
    T* Mgdy; //Mass matrix for cell with cell + dy
    T* Mgdxdy; //Mass matrix for cell with cell + dx + dy
    T* Mgmdxdy; //Mass matrix for cell with cell - dx + dy
    //Incomplete system matrix (without mass matrix terms): Nx*Ny*6*Nx*Ny*6
    //Complete system matrix (with mass matrix terms): Nx*Ny*6*Nx*Ny*6


    SpMat A; //Finite difference components of system matrix (6*Nx*Ny by 6*Nx*Ny)
    SpMat Am; //Mass matrix components of system matrix (6*Nx*Ny by 6*Nx*Ny)
    SpMat Ac; //Complete system matrix (6*Nx*Ny by 6*Nx*Ny)
    Eigen::VectorX<T> C; //System vector: Nx*Ny*6, same order as field

public:
    Field2D3V() = delete;
    explicit Field2D3V(const typename Config<T,2,3>::FieldConfig& fieldConfig,
                       const typename Config<T,2,3>::BCConfig& bcConfig);
    Field2D3V(const Field2D3V& other) = delete;
    ~Field2D3V();

    /*Initialisation works as following:
     * if initialiseFromSpecies is false, the field is only forcedE and forcedB values,
     * if initialiseFromSpecies is true then the field is solved to satisfy divergence laws and curl laws with all time
     * derivatives equal to zero, then forcedE and forcedB are superimposed (unless onlyForcedE or onlyForcedB are true
     * in which case they are only that). For now only works with fully periodic conditions.
     * It assumes that values in Z direction are zero for the derivatives, additionally, all average components of
     * electric field are 0 (before adding forced terms). */
    void initialise(const std::vector<Species<T,2,3>*>& species, T dt);

    //Save electric/magnetic field in 3 successive matrices, x, y and z components of the field. Rows are different x
    //values and columns are y values.
    void saveElectricField(std::ofstream& outputFile) const;
    void saveMagneticField(std::ofstream& outputFile) const;

    void saveEnergy(std::ofstream& outputFile) const;
    void saveElectrostaticPotential(std::ofstream &outputFile, const std::vector<Species<T,2,3>*> &species) const;
    void saveCurrent(std::ofstream& outputFile) const;

    const T* getField() const;
    const T* getFieldT() const;

    //Returns a pointer to a Nx by Ny matrix in col-major format (i.e. rows are different x-vals and cols are y-vals)
    //Electrostatic potential is defined to be 0 at 0,0
    std::unique_ptr<const T> getElectrostaticPotential(const std::vector<Species<T,2,3>*> &species) const;
private:
    void accumulateJ(const std::vector<Species<T,2,3>*>& species);
    void accumulateM(const std::vector<Species<T,2,3>*>& species, T dt);
    virtual void solveAndAdvance(T dt);
    virtual void initialiseSystemMatrix(T dt);
};


#endif //PIC_SEMI_IMPLICIT_FIELD2D3V_H