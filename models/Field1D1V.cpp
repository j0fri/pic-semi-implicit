#include <cmath>
#include <iostream>
#include "Field1D1V.h"

#define F77NAME(x) x##_
extern "C" {
void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B,
                    const int& ldb, int& info);
void F77NAME(sgesv)(const int& n, const int& nrhs, const float * A,
                    const int& lda, int * ipiv, float * B,
                    const int& ldb, int& info);
}


//TODO: tidy initialisation of nullptrs
template <typename T>
Field1D1V<T>::Field1D1V(const typename Config<T,1,1>::FieldConfig& fieldConfig): Field<T,1,1>(fieldConfig){
    Nx = fieldConfig.grid.dimensions[0].Nc;
    dx = fieldConfig.grid.getSpacings()[0];
}


//TODO: tidy initialisation of nullptrs and copy field values
template <typename T>
Field1D1V<T>::Field1D1V(const Field1D1V& other): Field<T,1,1>(other){
    Nx = other.Nx;
    dx = other.dx;
}

template <typename T>
Field1D1V<T>::~Field1D1V(){
    delete[] this->E;
    delete[] this->Et;
    delete[] this->J;
    delete[] this->Mgg;
    delete[] this->Mggp;
    delete[] this->A;
    delete[] this->C;
}

template<typename T>
void Field1D1V<T>::initialise(const std::vector<Species<T, 1, 1> *> &species) {
    this->E = new T[Nx];
    this->Et = new T[Nx];
    this->J = new T[Nx];
    this->Mgg = new T[Nx];
    this->Mggp = new T[Nx];
    this->A = new T[Nx*Nx];
    this->C = new T[Nx];

    T* rhoDist = new T[Nx];
    std::fill(rhoDist, rhoDist+Nx, 0);

    for(Species<T,1,1>* s: species){
        s->advancePositions((T)0, this);
        const int* g = s->getG();
        for(int i = 0; i < s->Np; ++i){
            rhoDist[g[i]] += s->q/dx;
        }
    }

    E[0] = 0;
    T shifting = 0;
    for(int i = 1; i < Nx; ++i){
        E[i] = E[i-1] + rhoDist[i]/this->e0*dx; //Not sure if *dx
        shifting += E[i]/Nx;
    }
    for(int i = 0; i < Nx; ++i){
        E[i] -= shifting;
    }

    delete[] rhoDist;

    std::fill(Et, Et+Nx, 0);
    std::fill(A, A+Nx*Nx, 0);
}

//TODO: clean and optimise
//TODO: ADD FORCED AND ADDED FIELD DISTRIBUTIONS
template <typename T>
void Field1D1V<T>::solveAndAdvance(T dt) {
    T c2 = 2/(this->c*dt);
    T c3 = 4*M_PI/this->c;
    for(unsigned int i = 0; i < Nx; ++i){
        A[i*(Nx+1)] = -c2;
        A[i*(Nx+1)] -= c3*Mgg[i];
        A[i*Nx+(i+1)%Nx] = -c3*Mggp[i];
        A[((i+1)%Nx)*Nx+i] = -c3*Mggp[i];
    }
    for(unsigned int i = 0; i < Nx; ++i){
        C[i] = -c2*E[i]+c3*J[i];
    }
    int* pivots = new int[Nx];
    int info = 0;

    if constexpr (std::is_same<T,float>::value){
        F77NAME(sgesv)((int)Nx, 1, A, (int)Nx, pivots, C, (int)Nx, info);
    }else if constexpr (std::is_same<T,double>::value){
        F77NAME(dgesv)((int)Nx, 1, A, (int)Nx, pivots, C, (int)Nx, info);
    }


    std::copy(C, C+Nx, Et);
    delete[] pivots;

    for(int i = 0; i < Nx; ++i){
        E[i] = 2*Et[i]-E[i];
    }
}

//TODO: clean
template<typename T>
void Field1D1V<T>::accumulateM(const std::vector<Species<T, 1, 1> *> &species, T dt) {
    //Clear Mgg and Mggp
    std::fill(this->Mgg, this->Mgg + Nx, 0);
    std::fill(this->Mggp, this->Mggp + Nx, 0);

    unsigned int Np;
    T m;
    T q;
    T beta;
    const int* g;
    const int* gp;
    const T* wg;
    const T* wgp;

    for(unsigned int s = 0; s < (unsigned int)species.size(); ++s){
        Np = species[s]->Np;
        m = species[s]->m;
        q = species[s]->q;
        beta = q/m*dt/2;
        g = species[s]->getG();
        gp = species[s]->getGp();
        wg = species[s]->getWg();
        wgp = species[s]->getWgp();
        for (unsigned int p = 0; p < Np; ++p){
            Mgg[g[p]] += q*beta*wg[p]*wg[p]/dx;
            Mgg[gp[p]] += q*beta*wgp[p]*wgp[p]/dx;
            Mggp[g[p]] += q*beta*wg[p]*wgp[p]/dx;
        }
    }
}

//TODO: clean
template<typename T>
void Field1D1V<T>::accumulateJ(const std::vector<Species<T, 1, 1> *> &species) {
    //Clear J 
    std::fill(this->J, this->J + Nx, 0);

    unsigned int Np;
    T q;
    const T* v;
    const int* g;
    const int* gp;
    const T* wg;
    const T* wgp;

    for (unsigned int s = 0; s < (unsigned int)species.size(); ++s){
        Np = species[s]->Np;
        q = species[s]->q;
        v = species[s]->getV();
        g = species[s]->getG();
        gp = species[s]->getGp();
        wg = species[s]->getWg();
        wgp = species[s]->getWgp();
        for (unsigned int p = 0; p < Np; ++p){
            J[g[p]] += q*v[p]*wg[p]/dx;
            J[gp[p]] += q*v[p]*wgp[p]/dx;
        }
    }
}

template<typename T>
const T *Field1D1V<T>::getEt() const {
    return static_cast<const T*>(Et);
}

template<typename T>
void Field1D1V<T>::saveElectricField(std::ofstream &outputFile) const {
    //TODO: implement this
}

template<typename T>
void Field1D1V<T>::saveMagneticField(std::ofstream &outputFile) const {
    //TODO: implement this
}

template<typename T>
void Field1D1V<T>::saveEnergy(std::ofstream &outputFile) const {
    T output = 0;
    for(unsigned int i = 0; i < Nx; ++i){
        output += E[i]*E[i];
    }
    outputFile << output/this->e0*dx/(8*M_PI) << std::endl;
}

template<typename T>
void Field1D1V<T>::saveVoltage(std::ofstream &outputFile) const {
    //TODO: implement this
}

template class Field1D1V<float>;
template class Field1D1V<double>;