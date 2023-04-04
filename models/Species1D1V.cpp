#include <algorithm>
#include <cmath>

#include "Species1D1V.h"
#include "Field1D1V.h"

template<typename T>
Species1D1V<T>::Species1D1V(const typename Config<T, 1, 1>::SpeciesConfig &speciesConfig,
                            const typename Config<T,1,1>::BCConfig& bcConfig):
    Species<T,1,1>(speciesConfig, bcConfig), g{nullptr}, gp{nullptr} {}

//TODO: handle null ptrs
template<typename T>
Species1D1V<T>::Species1D1V(const Species1D1V<T> &other): Species<T,1,1>(other) {
    //TODO: copy vectors and move constructor
}

template<typename T>
Species1D1V<T>::~Species1D1V() {
    delete[] this->x;
    delete[] this->v;
    delete[] this->g;
    delete[] this->gp;
    delete[] this->wg;
    delete[] this->wgp;
}

template <typename T>
const T* Species1D1V<T>::getV() const{
    return static_cast<const T*>(v);
}

template <typename T>
const int* Species1D1V<T>::getG() const{
    return static_cast<const int*>(g);
}

template <typename T>
const int* Species1D1V<T>::getGp() const{
    return static_cast<const int*>(gp);
}

template <typename T>
const T* Species1D1V<T>::getWg() const{
    return static_cast<const T*>(wg);
}

template <typename T>
const T* Species1D1V<T>::getWgp() const{
    return static_cast<const T*>(wgp);
}

template <typename T>
T Species1D1V<T>::getTotalKineticEnergy() const{
    float output = 0;
    for(int i = 0; i < this->Np; ++i){
        output += v[i]*v[i];
    }
    return output*this->m/2;
}

template <typename T>
void Species1D1V<T>::initialisePositions() {
    this->x = new T[this->Np];
    std::vector<std::array<T, 1>> tempX = this->initialXDist.generate(this->Np, this->initialXGrid);
    for(unsigned int i = 0; i < (unsigned int) tempX.size(); ++i){
        x[i] = tempX[i][0];
    }

    //TODO: override initialise to also call a method which does this:
    this->g = new int[this->Np];
    this->gp = new int[this->Np];
    this->wg = new T[this->Np];
    this->wgp = new T[this->Np];
}

template <typename T>
void Species1D1V<T>::initialiseVelocities(){
    this->v = new T[this->Np];
    std::vector<std::array<T, 1>> tempV = this->initialVDist.generate(this->Np, this->initialVGrid);
    for(unsigned int i = 0; i < (unsigned int) tempV.size(); ++i){
        v[i] = tempV[i][0];
    }
}

template <typename T>
void Species1D1V<T>::advancePositions(T dt, const Field<T, 1, 1>* field) {
    for(int i = 0; i < this->Np; ++i){
        x[i] += v[i]*dt;
        //Periodic boundary conditions
        while (x[i] >= field->grid.dimensions[0].max){
            x[i] -= field->grid.dimensions[0].max-field->grid.dimensions[0].min;
        }
        while (x[i] < field->grid.dimensions[0].min){
            x[i] += field->grid.dimensions[0].max-field->grid.dimensions[0].min;
        }
    }
    this->computeWeights(field);
    this->computeAlphas(field, dt);
}

template <typename T>
void Species1D1V<T>::advanceVelocities(T dt, const Field<T, 1, 1>* field) {
    const T* Et = ((Field1D1V<T>*)field)->getEt();
    T beta = this->q/this->m*dt/2;
    for(int i = 0; i < this->Np; ++i){
        v[i] += 2*beta*(Et[g[i]]*wg[i]+Et[gp[i]]*wgp[i]);
    }
}

//TODO: ADD THIS
template<typename T>
void Species1D1V<T>::initialisePositions(std::ifstream &file) {

}

//TODO: ADD THIS
template<typename T>
void Species1D1V<T>::initialiseVelocities(std::ifstream &file) {

}

template<typename T>
void Species1D1V<T>::computeAlphas(const Field<T, 1, 1> *field, T dt) {
    //TODO: tidy as it's not needed
}

template<typename T>
void Species1D1V<T>::computeWeights(const Field<T, 1, 1> *field) {
    T dx = field->grid.getSpacings()[0];
    unsigned int Nx = field->grid.dimensions[0].Nc;
    for(int i = 0; i < this->Np; ++i){
        g[i] = std::floor(x[i]/dx);
        gp[i] = (g[i]+1) % Nx;
        wg[i] = std::fmod(x[i], dx)/dx; //TODO: check it's not 1-this
        wgp[i] = 1-wg[i];
        if(g[i]<0 || g[i]>=Nx || gp[i]<0 || gp[i]>=Nx){
            throw std::runtime_error("Problem"); //TODO: clean or ignore
        }
    }
}

template<typename T>
void Species1D1V<T>::savePosition(std::ofstream &outputFile) const {
    for(unsigned int i = 0; i < this->Np; ++i){
        outputFile << x[i] << " ";
    }
    outputFile << std::endl;
}

template<typename T>
void Species1D1V<T>::savePositionDistribution(std::ofstream &outputFile, Field<T,1,1> *field) const {
    //TODO: implement this
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Species1D1V<T>::saveVelocity(std::ofstream &outputFile) const {
    for(unsigned int i = 0; i < this->Np; ++i){
        outputFile << v[i] << " ";
    }
    outputFile << std::endl;
}

template<typename T>
void Species1D1V<T>::saveVelocityDistribution(std::ofstream &outputFile) const {
    //TODO: implement this
}

template<typename T>
void Species1D1V<T>::saveEnergy(std::ofstream &outputFile) const {
    outputFile << this->getTotalKineticEnergy() << std::endl;
}


template class Species1D1V<float>;
template class Species1D1V<double>;