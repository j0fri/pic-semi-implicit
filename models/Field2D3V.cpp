//
// Created by jf1519 on 22/03/23.
//

#include "Field2D3V.h"

template<typename T>
Field2D3V<T>::Field2D3V(const typename Config<T, 2, 3>::FieldConfig &fieldConfig,
                        const typename Config<T,2,3>::BCConfig& bcConfig): Field<T,2,3>(fieldConfig, bcConfig),
               Ng{fieldConfig.grid.dimensions[0].Nc*fieldConfig.grid.dimensions[1].Nc},
               Nx{fieldConfig.grid.dimensions[0].Nc}, Ny{fieldConfig.grid.dimensions[1].Nc},
               field{new T[6*Ng]}, fieldT{new T[6*Ng]}, J{new T[3*Ng]},
               Mg{new T[9*Ng]}, Mgdx{new T[9*Ng]}, Mgdy{new T[9*Ng]}, Mgdxdy{new T[9*Ng]}, Mgmdxdy{new T[9*Ng]},
               A{new T[Ng*Ng*36]}, c{new T[Ng*6]} {}


template<typename T>
Field2D3V<T>::~Field2D3V() {
    delete[] field;
    delete[] fieldT;
    delete[] J;
    delete[] Mg;
    delete[] Mgdx;
    delete[] Mgdy;
    delete[] Mgdxdy;
    delete[] Mgmdxdy;
    delete[] A;
    delete[] c;
}

template<typename T>
void Field2D3V<T>::initialise(const std::vector<Species<T, 2, 3> *> &species) {
    if(!this->initialiseFromSpecies){
        std::fill(field, field+6*Ng,(T)0.0);
        T dx = this->grid.getSpacings()[0];
        T dy = this->grid.getSpacings()[1];
        T minx = this->grid.dimensions[0].min;
        T miny = this->grid.dimensions[1].max;
        for(unsigned int i = 0; i < Nx; ++i){
            for(unsigned int j = 0; j < Ny; ++j){
                for(unsigned int dim = 0; dim < 3; ++dim){
                    T tempE = this->forcedE[dim].f(std::array<T,2>{minx+i*dx,miny+j*dy});
                    T tempB = this->forcedB[dim].f(std::array<T,2>{minx+i*dx,miny+j*dy});
                    field[6*i+6*Nx*j+dim] = tempE;
                    field[6*i+6*Nx*j+dim+3] = tempB;
                }
            }
        }
    }else{
        //TODO: add divergence-free intialisation
        throw std::runtime_error("Not implemented yet.");
    }
}


template<typename T>
void Field2D3V<T>::saveElectricField(std::ofstream &outputFile) const {
    //TODO: add electric field save
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Field2D3V<T>::saveMagneticField(std::ofstream &outputFile) const {
    //TODO: add magnetic field save
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Field2D3V<T>::saveEnergy(std::ofstream &outputFile) const {
    //TODO: add field energy save
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Field2D3V<T>::saveVoltage(std::ofstream &outputFile) const {
    //TODO: add voltage save
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
const T *Field2D3V<T>::getField() const {
    return static_cast<const T*>(field);
}

template<typename T>
const T *Field2D3V<T>::getFieldT() const {
    return static_cast<const T*>(fieldT);
}


template class Field2D3V<float>;
template class Field2D3V<double>;