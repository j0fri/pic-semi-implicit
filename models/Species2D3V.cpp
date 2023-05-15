#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Species2D3V.h"
#include "Field2D3V.h"
#include "../helpers/math_helper.h"
#include "../helpers/output_helper.h"

template<typename T>
Species2D3V<T>::Species2D3V(const typename Config<T,2,3>::SpeciesConfig &speciesConfig,
                            const  typename Config<T,2,3>::BCConfig& bcConfig)
                            :Species<T,2,3>(speciesConfig, bcConfig),
                            pos{Vector2<T>(this->Np)}, vel{Vector3<T>(this->Np)},
                            g{Vector2<unsigned int>(this->Np)}, gp{Vector2<unsigned int>(this->Np)},
                            gB{Vector2<unsigned int>(this->Np)}, gpB{Vector2<unsigned int>(this->Np)},
                            wg{Vector2<T>(this->Np)}, wgp{Vector2<T>(this->Np)},
                            wgB{Vector2<T>(this->Np)}, wgpB{Vector2<T>(this->Np)},
                            Ep{Vector3<T>(this->Np)}, Bp{Vector3<T>(this->Np)},
                            alpha{nullptr}, vBar{Vector3<T>(this->Np)}{
    try{
        alpha = new T[9*this->Np];
    }catch(const std::bad_alloc& e){
        std::cerr << "Exception thrown in species memory allocation O(Np)" << std::endl;
        throw;
    }
}

template<typename T>
Species2D3V<T>::Species2D3V(const Species2D3V<T> &other)=default;

template<typename T>
Species2D3V<T>::~Species2D3V(){
    delete[] alpha;
};

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

template<typename T>
const T *Species2D3V<T>::getAlpha() const {
    return static_cast<const T*>(alpha);
}

template<typename T>
T Species2D3V<T>::getTotalKineticEnergy() const {
    T total = 0;
    for(unsigned int i = 0; i < 3*this->Np; ++i){
        total += vel.x[i]*vel.x[i];
    }
    T kineticEnergy = total/2*this->m;
    return kineticEnergy;
}

template<typename T>
void Species2D3V<T>::advancePositions(T dt, const Field<T,2,3> *field) {
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
            while(ptr[i]>=max){
                ptr[i] -= length;
            }
            while(ptr[i]<min){
                ptr[i] += length;
            }
        }
    }
    this->computeWeights(field);
    this->computeLocalB(field);
    this->computeAlphas(dt);
}

template<typename T>
void Species2D3V<T>::advanceVelocities(T dt, const Field<T, 2, 3> *field) {
    this->computeLocalE(field);
    Vector3<T> temp = vel + Ep*(this->q*dt/(2*this->m));
    for(unsigned int i = 0; i < this->Np; ++i){
        math_helper::gemv(3,3,(T)1.0,alpha+9*i,3,temp.x+i,temp.n,vBar.x+i,vBar.n);
    }
    vel = vBar*((T)2.0) - vel;
}

template<typename T>
void Species2D3V<T>::savePosition(std::ofstream &outputFile) const {
    if(this->Np > 1000){
        std::cout << "WARNING: consider disabling savePosition if Np > 1000" << std::endl;
    }
    //outputFile << this->m << " " << this->q << std::endl;
    output_helper::outputRowMajorMatrix(this->pos.x,2,this->Np,this->Np,1,outputFile);
}

template<typename T>
void Species2D3V<T>::savePositionDistribution(std::ofstream &outputFile, Field<T,2,3> *field) const {
    outputFile << this->m << " " << this->q << std::endl;
    unsigned int Nx = field->grid.dimensions[0].Nc;
    unsigned int Ny = field->grid.dimensions[1].Nc;
    auto freq = new unsigned int[Nx*Ny];
    std::fill(freq,freq+Nx*Ny,0);
    for(unsigned int p = 0; p < this->Np; ++p){
        ++freq[g.x[p]+g.y[p]*Nx];
    }
    output_helper::outputColMajorMatrix(freq,(int)Nx,(int)Ny,(int)Nx,1,outputFile);
    delete[] freq;
}

template<typename T>
void Species2D3V<T>::saveVelocity(std::ofstream &outputFile) const {
    if(this->Np > 1000){
        std::cout << "WARNING: consider disabling saveVelocity if Np > 1000" << std::endl;
    }
    //outputFile << this->m << " " << this->q << std::endl;
    output_helper::outputRowMajorMatrix(this->vel.x,3,this->Np,this->Np,1,outputFile);
}

template<typename T>
void Species2D3V<T>::saveVelocityDistribution(std::ofstream &outputFile) const {
    //TODO: implement
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Species2D3V<T>::saveEnergy(std::ofstream &outputFile) const {
    outputFile << this->m << " " << this->q << std::endl;
    outputFile << this->getTotalKineticEnergy() << std::endl;
}

template<typename T>
void Species2D3V<T>::initialisePositions() {
    std::vector<std::array<T, 2>> tempPos = this->initialXDist.generate(this->Np, this->initialXGrid);
    for(unsigned int i = 0; i < (unsigned int) tempPos.size(); ++i){
        pos.x[i] = tempPos[i][0];
        pos.y[i] = tempPos[i][1];
    }
}

template<typename T>
void Species2D3V<T>::initialiseVelocities() {
    std::vector<std::array<T, 3>> tempVel = this->initialVDist.generate(this->Np, this->initialVGrid);
    for(unsigned int i = 0; i < (unsigned int) tempVel.size(); ++i){
        vel.x[i] = tempVel[i][0];
        vel.y[i] = tempVel[i][1];
        vel.z[i] = tempVel[i][2];
    }
}

template<typename T>
void Species2D3V<T>::initialisePositions(std::ifstream &file) {
    try{
        if(!file.is_open()){
            throw std::invalid_argument("Initial position file not open.");
        }
        for(auto ptr = pos.x; ptr < pos.x + 2*this->Np; ++ptr){
            if(file.eof()){
                throw std::invalid_argument("Initial position file contains less columns than Np.");
            }
            file >> *ptr;
        }
    }catch(const std::exception& e){
        std::cout << "Initial position file has incorrect format." << std::endl;
        throw;
    }
}

template<typename T>
void Species2D3V<T>::initialiseVelocities(std::ifstream &file) {
    try{
        if(!file.is_open()){
            throw std::invalid_argument("Initial velocity file not open.");
        }
        for(auto ptr = vel.x; ptr < vel.x + 3*this->Np; ++ptr){
            if(file.eof()){
                throw std::invalid_argument("Initial velocity file contains less columns than Np.");
            }
            file >> *ptr;
        }
    }catch(const std::exception& e){
        std::cout << "Initial velocity file has incorrect format." << std::endl;
        throw;
    }
}

template<typename T>
void Species2D3V<T>::computeWeights(const Field<T,2,3>* field) {
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
        T invspacing = (T)1/spacing;
        T halfspacing = spacing/2;
        T shift1 = max-2*min-halfspacing;
        T shift2 = max-min;
        T shift3 = -min-halfspacing;
        unsigned int temp;
        T temp2;
        //Periodic base case:
        for(unsigned int i = 0; i < this->Np; ++i){
            //Unoptimised code:
//            gPtr[i] = (unsigned int) std::floor((posPtr[i]-min) / spacing);
//            gpPtr[i] = (gPtr[i]+1) % Nc;
//            gBPtr[i] = (((unsigned int) std::floor((posPtr[i]+max-2*min-spacing/2) / spacing))) % Nc;
//            gpBPtr[i] = (gBPtr[i]+1) % Nc;
//            wgpPtr[i] = std::fmod(posPtr[i]+max-2*min, spacing)/spacing;
//            wgPtr[i] = (T)1-wgpPtr[i];
//            wgpBPtr[i] = std::fmod(posPtr[i]+max-2*min-spacing/2, spacing)/spacing;
//            wgBPtr[i] = (T)1-wgpBPtr[i];
//            std::cout << "pos: " << posPtr[i] << ", g: " << gPtr[i] << ", gp: " << gpPtr[i] << ", gB: " << gBPtr[i] << ", gpB: " << gpBPtr[i] << std::endl;


            gPtr[i] = (unsigned int)((posPtr[i]-min) * invspacing);
            temp = gPtr[i] + 1;
            gpPtr[i] = temp == Nc ? 0 : temp;
            temp = (unsigned int) ((posPtr[i]+shift1) * invspacing);
            gBPtr[i] = temp > Nc ? temp - Nc : temp;
            gpBPtr[i] = gBPtr[i] == Nc-1 ? 0 : gBPtr[i]+1;

            wgpPtr[i] = (posPtr[i]-(min+gPtr[i]*spacing)) * invspacing;
            wgPtr[i] = (T)1-wgpPtr[i];
            temp2 = posPtr[i]+shift3;
            temp2 = temp2 < 0 ? temp2 + shift2 : temp2;
            wgpBPtr[i] = (temp2-gBPtr[i]*spacing) * invspacing;
            wgBPtr[i] = (T)1-wgpBPtr[i];

        }
        if(!field->bcConfig.periodic[dim]){ //Fix non-periodicity //TODO: TEST
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
    unsigned int i1, i2, i3, i4;
    T w1, w2, w3, w4;
    unsigned int Nx = field->grid.dimensions[0].Nc;
    for(unsigned int i = 0; i < this->Np; ++i){
        i1 = gB.x[i] + Nx * gB.y[i];
        i2 = gpB.x[i] + Nx * gB.y[i];
        i3 = gB.x[i] + Nx * gpB.y[i];
        i4 = gpB.x[i] + Nx * gpB.y[i];
        w1 = wgB.x[i]*wgB.y[i];
        w2 = wgpB.x[i]*wgB.y[i];
        w3 = wgB.x[i]*wgpB.y[i];
        w4 = wgpB.x[i]*wgpB.y[i];
        for(unsigned int dim = 3; dim < 6; ++dim){
            Bp[dim-3][i] = w1*f[6*i1+dim] + w2*f[6*i2+dim] + w3*f[6*i3+dim] + w4*f[6*i4+dim];
        }
    }
}

template<typename T>
void Species2D3V<T>::computeAlphas(T dt) {
    //From Giovanni Lapenta:
    //    function alpha = alpha(beta,Bx,By,Bz )
    //                     %B=[Bx By Bz];
    //    %I=diag([1 1 1]);
    //    %IcB=[cross(I(1,:),B); cross(I(2,:),B); cross(I(3,:),B)];
    //
    //    %alpha=(I-beta*IcB+beta^2*B'*B)/(1+beta^2*B*B');
    //
    //    sx=Bx*beta;sy=By*beta;sz=Bz*beta;
    //    alpha=[1+sx*sx  sz+sx*sy   -sy+sx*sz;
    //    -sz+sx*sy  1+sy*sy   sx+sy*sz;
    //    sy+sx*sz   -sx+sy*sz    1+sz*sz]/(1+sx*sx+sy*sy+sz*sz);
    //    end
    unsigned int ai;
    T beta = dt*this->q / (2*this->m);
    for(unsigned int i = 0; i < this->Np; ++i){
        ai = 9*i;
        T sx = Bp.x[i] * beta;
        T sy = Bp.y[i] * beta;
        T sz = Bp.z[i] * beta;
        alpha[ai]=1+sx*sx;
        alpha[ai+1]=-sz+sx*sy;
        alpha[ai+2]=sy+sx*sz;
        alpha[ai+3]=sz+sx*sy;
        alpha[ai+4]=1+sy*sy;
        alpha[ai+5]=-sx+sy*sz;
        alpha[ai+6]=-sy+sx*sz;
        alpha[ai+7]=sx+sy*sz;
        alpha[ai+8]=1+sz*sz;
        T norm = (T)1/(1+sx*sx+sy*sy+sz*sz);
        for(unsigned int j = 0; j < 9; ++j){
            alpha[ai+j]*=norm;
        }
    }
}

template<typename T>
void Species2D3V<T>::computeLocalE(const Field<T, 2, 3> *field) {
    const T* f = ((Field2D3V<T>*)field)->getFieldT();
    unsigned int i1, i2, i3, i4;
    T w1, w2, w3, w4;
    unsigned int Nx = field->grid.dimensions[0].Nc;
    for(unsigned int i = 0; i < this->Np; ++i){
        i1 = g.x[i] + Nx*g.y[i];
        i2 = gp.x[i] + Nx*g.y[i];
        i3 = g.x[i] + Nx*gp.y[i];
        i4 = gp.x[i] + Nx*gp.y[i];
        w1 = wg.x[i]*wg.y[i];
        w2 = wgp.x[i]*wg.y[i];
        w3 = wg.x[i]*wgp.y[i];
        w4 = wgp.x[i]*wgp.y[i];
        for(unsigned int dim = 0; dim < 3; ++dim){
            Ep[dim][i] = w1*f[6*i1+dim] + w2*f[6*i2+dim] + w3*f[6*i3+dim] + w4*f[6*i4+dim];
        }
    }
}


template class Species2D3V<float>;
template class Species2D3V<double>;