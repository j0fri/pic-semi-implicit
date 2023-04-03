#include <cmath>
#include <iostream>
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

//TODO: ADD KINETIC ENERGY
template<typename T>
T Species2D3V<T>::getTotalKineticEnergy() const {
    return 0;
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
//    std::ofstream file1("outputs/vBar.txt", std::ios::app);
//    file1 << "Begin advanceVelocities: " << std::endl;
//    output_helper::outputRowMajorMatrix(vel.x, 3, this->Np,this->Np,1,file1);


    this->computeLocalE(field);


//    auto temp2 = (Ep*(this->q*dt/(2*this->m)));
//    file1 << "Acceleration term: " << std::endl;
//    output_helper::outputRowMajorMatrix(temp2.x, 3, this->Np,this->Np,1,file1);


    Vector3<T> temp = vel + Ep*(this->q*dt/(2*this->m));
    for(unsigned int i = 0; i < this->Np; ++i){
        math_helper::gemv(3,3,(T)1.0,alpha+9*i,3,temp.x+i,temp.n,vBar.x+i,vBar.n);
    }
    vel = vBar*((T)2.0) - vel;

//    file1 << "Local magnetic field: " << std::endl;
//    output_helper::outputRowMajorMatrix(Bp.x,3,this->Np,this->Np,1,file1);
//    file1 << "alpha1: " << std::endl;
//    output_helper::outputColMajorMatrix(alpha,3,3,3,1,file1);
//    file1 << "alpha2: " << std::endl;
//    output_helper::outputColMajorMatrix(alpha+9,3,3,3,1,file1);
//    file1 << "alpha3: " << std::endl;
//    output_helper::outputColMajorMatrix(alpha+18,3,3,3,1,file1);
//    file1 << "alpha4: " << std::endl;
//    output_helper::outputColMajorMatrix(alpha+27,3,3,3,1,file1);
//    file1 << "alpha5: " << std::endl;
//    output_helper::outputColMajorMatrix(alpha+36,3,3,3,1,file1);
//    file1 << "vbar: " << std::endl;
//    output_helper::outputRowMajorMatrix(vBar.x, 3, this->Np,this->Np,1,file1);
//    file1 << "Finish advanceVelocities: " << std::endl;
//    output_helper::outputRowMajorMatrix(vel.x, 3, this->Np,this->Np,1,file1);
//    file1 << std::endl << std::endl << std::endl;
}

template<typename T>
void Species2D3V<T>::savePosition(std::ofstream &outputFile) const {
    if(this->Np > 1000){
        std::cout << "WARNING: consider disabling savePosition if Np > 1000" << std::endl;
    }
    outputFile << this->m << " " << this->q << std::endl;
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
        std::cout << "WARNING: consider disabling savePosition if Np > 1000" << std::endl;
    }
    outputFile << this->m << " " << this->q << std::endl;
    output_helper::outputRowMajorMatrix(this->vel.x,3,this->Np,this->Np,1,outputFile);
}

template<typename T>
void Species2D3V<T>::saveVelocityDistribution(std::ofstream &outputFile) const {
    //TODO: implement
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Species2D3V<T>::saveEnergy(std::ofstream &outputFile) const {
    //TODO: implement
    throw std::runtime_error("Not implemented yet.");
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
void Species2D3V<T>::initialisePositions(const std::ifstream &file) {
    //TODO: initialisation from file
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Species2D3V<T>::initialiseVelocities(const std::ifstream &file) {
    //TODO: initialisation from file
    throw std::runtime_error("Not implemented yet.");
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
        //Periodic base case:
        for(unsigned int i = 0; i < this->Np; ++i){
            gPtr[i] = (unsigned int) std::floor((posPtr[i]-min) / spacing);
            gpPtr[i] = (gPtr[i]+1) % Nc;
            gBPtr[i] = (((unsigned int) std::floor((posPtr[i]-min-spacing/2) / spacing)) + Nc) % Nc;
            gpBPtr[i] = (gBPtr[i]+1) % Nc;
            wgpPtr[i] = std::fmod(posPtr[i]+max, spacing)/spacing;
            wgPtr[i] = (T)1-wgpPtr[i];
            wgpBPtr[i] = std::fmod(posPtr[i]+max-spacing/2, spacing)/spacing;
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
    for(unsigned int i = 0; i < this->Np; ++i){
        unsigned int i1 = gB.x[i] + field->grid.dimensions[0].Nc * gB.y[i];
        unsigned int i2 = gpB.x[i] + field->grid.dimensions[0].Nc * gB.y[i];
        unsigned int i3 = gB.x[i] + field->grid.dimensions[0].Nc * gpB.y[i];
        unsigned int i4 = gpB.x[i] + field->grid.dimensions[0].Nc * gpB.y[i];
        T w1 = wgB.x[i]*wgB.y[i];
        T w2 = wgpB.x[i]*wgB.y[i];
        T w3 = wgB.x[i]*wgpB.y[i];
        T w4 = wgpB.x[i]*wgpB.y[i];
        for(unsigned int dim = 0; dim < 3; ++dim){
            Bp[dim][i] = w1*f[6*i1+3+dim] + w2*f[6*i2+3+dim] + w3*f[6*i3+3+dim] + w4*f[6*i4+3+dim];
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
            Ep[dim][i] = w1*f[6*i1+dim] + w2*f[6*i2+dim] + w3*f[6*i3+dim] + w4*f[6*i4+dim];
        }
    }
}


template class Species2D3V<float>;
template class Species2D3V<double>;