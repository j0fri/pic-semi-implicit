//
// Created by jf1519 on 22/03/23.
//

#include <cmath>
#include "Species2D3V.h"
#include "Field2D3V.h"

template<typename T>
Species2D3V<T>::Species2D3V(const Config<T,2,3>::SpeciesConfig &speciesConfig, const Config<T,2,3>::BCConfig& bcConfig)
        :Species<T,2,3>(speciesConfig, bcConfig),
        pos{Vector2<T>(this->Np)}, vel{Vector3<T>(this->Np)},
        g{Vector2<unsigned int>(this->Np)}, gp{Vector2<unsigned int>(this->Np)},
        gB{Vector2<unsigned int>(this->Np)}, gpB{Vector2<unsigned int>(this->Np)},
        wg{Vector2<T>(this->Np)}, wgp{Vector2<T>(this->Np)},
        wgB{Vector2<T>(this->Np)}, wgpB{Vector2<T>(this->Np)},
        Ep{Vector3<T>(this->Np)}, Bp{Vector3<T>(this->Np)}{
}

template<typename T>
Species2D3V<T>::Species2D3V(const Species2D3V<T> &other)=default;

template<typename T>
Species2D3V<T>::~Species2D3V()=default;

template<typename T>
const Vector3<T>& Species2D3V<T>::getVel() const {
    return vel;
}

template<typename T>
const Vector2<unsigned int>& Species2D3V<T>::getG() const {
    return g;
}

template<typename T>
const Vector2<unsigned int>& Species2D3V<T>::getGp() const {
    return gp;
}

template<typename T>
const Vector2<T>& Species2D3V<T>::getWg() const {
    return wg;
}

template<typename T>
const Vector2<T>& Species2D3V<T>::getWgp() const {
    return wgp;
}

//TODO: ADD KINETIC ENERGY
template<typename T>
T Species2D3V<T>::getTotalKineticEnergy() const {
    return 0;
}

template<typename T>
void Species2D3V<T>::advancePositions(T dt, const Field<T, 2, 3> *field) {
    //Algorithm step 1
    //TODO: non-periodic boundary conditions + optimise
    pos += (vel * dt); //Advance vector
    //Periodic boundary conditions:
    for(unsigned int dim = 0; dim < 2; ++dim){
        if(!this->bcConfig.periodic[dim]){
            throw std::runtime_error("Non-boundary conditions not implemented yet.");
        }
        T* ptr = dim == 0 ? pos.x : pos.y;
        T min = field->grid.dimensions[dim].min;
        T max = field->grid.dimensions[dim].max;
        T length = max - min;
        for(unsigned int i = 0; i < this->Np; ++i){
            while(ptr[i]>max){
                ptr[i] -= length;
            }
            while(ptr[i]<min){
                ptr[i] += length;
            }
        }
    }
    this->computeWeights(field);
    this->computeLocalB(field);
    this->computeAlphas(field);
}

template<typename T>
void Species2D3V<T>::computeWeights(const Field<T, 2, 3> *field) {
    for(unsigned int dim = 0; dim < 2; ++dim){
        T min = field->grid.dimensions[dim].min;
        T max = field->grid.dimensions[dim].max;
        unsigned int Nc = field->grid.dimensions[dim].Nc;
        T spacing = field->grid.getSpacings()[dim];
        unsigned int* gPtr = g[dim];
        unsigned int* gpPtr = gp[dim];
        unsigned int* gBPtr = gB[dim];
        unsigned int* gpBPtr = gpB[dim];
        T* wgPtr = wg[dim];
        T* wgpPtr = wgp[dim];
        T* wgBPtr = wgB[dim];
        T* wgpBPtr = wgpB[dim];
        T* posPtr = pos[dim];
        //Periodic base case:
        for(unsigned int i = 0; i < this->Np; ++i){
            gPtr[i] = (unsigned int) std::floor((posPtr[i]-min) / spacing);
            gpPtr[i] = (gPtr[i]+1) % Nc;
            gBPtr[i] = (((unsigned int) std::floor((posPtr[i]-min-spacing/2) / spacing)) + Nc) % Nc;
            gpBPtr[i] = (gBPtr[i]+1) % Nc;
            wgpPtr[i] = std::fmod(posPtr[i], spacing)/spacing;
            wgPtr[i] = 1.0-wgpPtr[i];
            wgpBPtr[i] = std::fmod(posPtr[i]-spacing/2, spacing)/spacing;
            wgPtr[i] = 1.0-wgBPtr[i];
        }
        if(!field->bcConfig.periodic[dim]){ //Fix non-periodicity
            for(unsigned int i = 0; i < this->Np; ++i){
                if(posPtr[i]-min < spacing/2){
                    wgBPtr[i] = 0;
                }
                if(max-posPtr[i] <= spacing){
                    wgpPtr[i] = 0;
                    if(max-posPtr[i] <= spacing/2){
                        wgpBPtr[i] = 0;
                    }
                }
            }
        }
    }
}

template<typename T>
void Species2D3V<T>::computeLocalB(const Field<T, 2, 3> *field) {
    const T* f = ((Field2D3V<T>*)field)->getField();
    for(unsigned int i = 0; i < this->Np; ++i){
        unsigned int i1 = g.x[i] + field->grid.dimensions[0].Nc * g.y[i];
        unsigned int i2 = gp.x[i] + field->grid.dimensions[0].Nc * g.y[i];
        unsigned int i3 = g.x[i] + field->grid.dimensions[0].Nc * gp.y[i];
        unsigned int i4 = gp.x[i] + field->grid.dimensions[0].Nc * gp.y[i];
        T w1 = wg.x[i]*wg.y[i];
        T w2 = wgp.x[i]*wg.y[i];
        T w3 = wg.x[i]*wgp.y[i];
        T w4 = wgp.x[i]*wgp.y[i];
        for(unsigned int dim = 0; dim < 3; ++dim){
            Bp[dim][i] = w1*f[6*i1+3+dim] + w2*f[6*i2+3+dim] + w3*f[6*i3+3+dim] + w4*f[6*i4+3+dim];
        }
    }
}

template class Species2D3V<float>;
template class Species2D3V<double>;
