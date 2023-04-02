#include <iostream>

#include "Field2D3V.h"
#include "Species2D3V.h"
#include "../helpers/math_helper.h"
#include "../helpers/output_helper.h"

#define F77NAME(x) x##_
extern "C" {
void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B,
                    const int& ldb, int& info);
void F77NAME(sgesv)(const int& n, const int& nrhs, const float * A,
                    const int& lda, int * ipiv, float * B,
                    const int& ldb, int& info);
}


template<typename T>
Field2D3V<T>::Field2D3V(const typename Config<T,2,3>::FieldConfig &fieldConfig,
                        const typename Config<T,2,3>::BCConfig& bcConfig):
                        Field<T,2,3>(fieldConfig, bcConfig),
                        Ng{fieldConfig.grid.dimensions[0].Nc*fieldConfig.grid.dimensions[1].Nc},
                        Nx{fieldConfig.grid.dimensions[0].Nc}, Ny{fieldConfig.grid.dimensions[1].Nc}{
    try{
        field = new T[6*Ng];
        fieldT = new T[6*Ng];
        J = new T[3*Ng];
        Mg = new T[9*Ng];
        Mgdx = new T[9*Ng];
        Mgdy = new T[9*Ng];
        Mgdxdy = new T[9*Ng];
        Mgmdxdy = new T[9*Ng];
    }catch(const std::bad_alloc& e){
        std::cerr << "Error during field memory allocation of O(Nx*Ny) variables." << std::endl;
        throw;
    }

    try{
        A = new T[36*Ng*Ng];
        c = new T[6*Ng];
    }catch(const std::bad_alloc& e){
        std::cerr << "Error during field memory allocation of O(Nx*Ny*Nx*Ny) variables." << std::endl;
        throw;
    }
}


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
    }else if(this->bcConfig.periodic[0] && this->bcConfig.periodic[1]){ //Periodic in both dimensions
        T dx = this->grid.getSpacings()[0];
        T dy = this->grid.getSpacings()[1];

        std::fill(field, field+6*Ng,(T)0.0);

        std::unique_ptr<const T> phi = this->getElectrostaticPotential(species);
        T minx = this->grid.dimensions[0].min;
        T miny = this->grid.dimensions[1].max;

        unsigned int gdxi, gmdxi, gdyi, gmdyi;
        for(unsigned int i = 0; i < Nx; ++i){
            for(unsigned int j = 0; j < Ny; ++j){
                //Forced electric terms
                for(unsigned int dim = 0; dim < 3; ++dim){
                    T tempE = this->forcedE[dim].f(std::array<T,2>{minx+i*dx,miny+j*dy});
                    field[6*i+6*Nx*j+dim] = tempE;
                }
                //Forced magnetic terms
                for(unsigned int dim = 0; dim < 3; ++dim){
                    T tempB = this->forcedB[dim].f(std::array<T,2>{minx+i*dx,miny+j*dy});
                    field[6*i+6*Nx*j+3+dim] = tempB;
                }
                //Species electric terms
                gdxi = j*Nx+((i+1)%Nx);
                gmdxi = j*Nx+((i-1+Nx)%Nx);
                gdyi = ((j+1)%Ny)*Nx+i;
                gmdyi = ((j-1+Ny)%Ny)*Nx+i;
                if(!this->onlyForcedE){
                    field[6*i+6*Nx*j] -= (phi.get()[gdxi]-phi.get()[gmdxi])/(2*dx); //Ex = -dPhi/dx
                    field[6*i+6*Nx*j+1] -= (phi.get()[gdyi]-phi.get()[gmdyi])/(2*dy); //Ey = -dPhi/dy
                }
            }
        }

        //TODO: add divergence-free initialisation for magnetic field
    }else{
        //TODO: add divergence-free initialisation for non-periodic boundary conditions
        throw std::invalid_argument("Divergence-free initialisation for non-periodic boundary conditions not implemented"
                                    " yet.");
    }
}


template<typename T>
void Field2D3V<T>::saveElectricField(std::ofstream &outputFile) const {
    output_helper::outputColMajorMatrix(field,Nx,Ny,Nx*6,6,outputFile);
    output_helper::outputColMajorMatrix(field+1,Nx,Ny,Nx*6,6,outputFile);
    output_helper::outputColMajorMatrix(field+2,Nx,Ny,Nx*6,6,outputFile);
}

template<typename T>
void Field2D3V<T>::saveMagneticField(std::ofstream &outputFile) const {
    output_helper::outputColMajorMatrix(field+3,Nx,Ny,Nx*6,6,outputFile);
    output_helper::outputColMajorMatrix(field+4,Nx,Ny,Nx*6,6,outputFile);
    output_helper::outputColMajorMatrix(field+5,Nx,Ny,Nx*6,6,outputFile);
}

template<typename T>
void Field2D3V<T>::saveEnergy(std::ofstream &outputFile) const {
    //TODO: add field energy save
    throw std::runtime_error("Not implemented yet.");
}

template<typename T>
void Field2D3V<T>::saveElectrostaticPotential(std::ofstream &outputFile, const std::vector<Species<T,2,3>*> &species) const {
    std::unique_ptr<const T> phi = this->getElectrostaticPotential(species);
    output_helper::outputColMajorMatrix(&*phi,Nx,Ny,Nx,1,outputFile);
}

template<typename T>
const T *Field2D3V<T>::getField() const {
    return static_cast<const T*>(field);
}

template<typename T>
const T *Field2D3V<T>::getFieldT() const {
    return static_cast<const T*>(fieldT);
}

template<typename T>
std::unique_ptr<const T> Field2D3V<T>::getElectrostaticPotential(const std::vector<Species<T,2,3>*> &species) const {
    if(!(this->bcConfig.periodic[0] && this->bcConfig.periodic[1])){
        //TODO: add electrostatic potential for non-periodic boundary conditions
        throw std::runtime_error("Electrostatic potential not implemented for non-periodic boundary conditions");
    }

    T dx = this->grid.getSpacings()[0];
    T dy = this->grid.getSpacings()[1];

    T* EA;
    T* Ec;
    T* Q;
    int* pivots;

    try{
        EA = new T[Ng*Ng]; //System matrix for phi
        Ec = new T[Ng]; //System vector for phi, later becomes electrostatic potential at electric field grid
        Q = new T[Ng]; //Vector of charge in every cell at magnetic field grid
        pivots = new int[Ng]; //Pivot vector for matrix solve
    }catch(const std::bad_alloc& e){
        std::cerr << "Error in field memory temporary allocation in initialisation (O(Nx*Ny*Nx*Ny))" << std::endl;
        throw;
    }

    //Accumulate charges:
    std::fill(Q,Q+this->Ng,0.0);
    for(auto sPtr: species){
        const Vector2<unsigned int> gB = ((Species2D3V<T>*)sPtr)->getG();
        T q = sPtr->q;
        unsigned int gi;
        for(unsigned int i = 0; i < sPtr->Np; ++i){
            gi = gB.x[i]+this->Nx*gB.y[i];
            Q[gi] += q;
        }
    }

    //Construct system vector:
    for(unsigned int g = 0; g < this->Ng; ++g){
        Ec[g] = -Q[(g-1-Nx+Ng)%Ng]/(dx*dy*this->e0);
    }
    delete[] Q;

    //Populate system matrix:
    std::fill(EA, EA+Ng*Ng,0.0);
    unsigned int eq;
    unsigned int gi, gdxi, gmdxi, gdyi, gmdyi;
    T c1 = (T)1/(dx*dx);
    T c2 = (T)1/(dy*dy);
    for(unsigned int j = 0; j < this->Ny; ++j){
        for(unsigned int i = 0; i < this->Nx; ++i){
            eq = j*this->Nx+i;
            if(eq == 0){
                continue; //This will be the fixed point
            }
            gi = j*Nx+i;
            gdxi = j*Nx+((i+1)%Nx);
            gmdxi = j*Nx+((i-1+Nx)%Nx);
            gdyi = ((j+1)%Ny)*Nx+i;
            gmdyi = ((j-1+Ny)%Ny)*Nx+i;

            //Poisson equation:
            EA[gdxi*Ng+eq] += c1;
            EA[gmdxi*Ng+eq] += c1;
            EA[gdyi*Ng+eq] += c2;
            EA[gmdyi*Ng+eq] += c2;
            EA[gi*Ng+eq] += -2*c1-2*c2;
        }
    }

    //Potential is 0 at 0,0 to avoid matrix rank-deficiency
    EA[0] = 1;
    Ec[0] = (T)0;

    //Solve matrix
    if constexpr(std::is_same<T,double>::value){
        int info;
        F77NAME(dgesv)(this->Ng,1,(double *)EA,this->Ng,pivots, (double *)Ec, this->Ng,info);
        if(info != 0){
            throw std::runtime_error("System solve error during initialisation of electric field.");
        }
        delete[] pivots;
    }else{
        //TODO: Make non-lapack dependent
        throw std::invalid_argument("T can only be double when using Lapack for initialisation");
    }

    delete[] EA;

    return std::unique_ptr<const T>(Ec);
}

template<typename T>
void Field2D3V<T>::accumulateJ(const std::vector<Species<T,2,3>*> &species) {
    std::fill(J,J+3*Ng,(T)0.0);
    for(auto sPtr: species){
        T dJ[3];
        const Vector2<unsigned int> g = ((Species2D3V<T>*)sPtr)->getG();
        const Vector2<unsigned int> gp = ((Species2D3V<T>*)sPtr)->getGp();
        const Vector2<T> wg = ((Species2D3V<T>*)sPtr)->getWg();
        const Vector2<T> wgp = ((Species2D3V<T>*)sPtr)->getWgp();
        const Vector3<T> vel = ((Species2D3V<T>*)sPtr)->getVel();
        const T* alpha = ((Species2D3V<T>*)sPtr)->getAlpha();
        for(unsigned int p = 0; p < sPtr->Np; ++p){
            //Bottom-left cell
            unsigned int gi = g.x[p]+g.y[p]*Nx;
            math_helper::gemv(3,3,(T)(sPtr->q*wg.x[p]*wg.y[p]),alpha+9*gi,3,vel.x+p,sPtr->Np,dJ,1);
            for(unsigned int dim = 0; dim < 3; ++dim){
                J[3*gi+dim] +=dJ[dim];
            }
            //Bottom-right cell (dx)
            gi = gp.x[p]+g.y[p]*Nx;
            math_helper::gemv(3,3,(T)(sPtr->q*wgp.x[p]*wg.y[p]),alpha+9*gi,3,vel.x+p,sPtr->Np,dJ,1);
            for(unsigned int dim = 0; dim < 3; ++dim){
                J[3*gi+dim] +=dJ[dim];
            }
            //Top-left cell (dy)
            gi = g.x[p]+gp.y[p]*Nx;
            math_helper::gemv(3,3,(T)(sPtr->q*wg.x[p]*wgp.y[p]),alpha+9*gi,3,vel.x+p,sPtr->Np,dJ,1);
            for(unsigned int dim = 0; dim < 3; ++dim){
                J[3*gi+dim] +=dJ[dim];
            }
            //Top-right cell (dx, dy)
            gi = gp.x[p]+gp.y[p]*Nx;
            math_helper::gemv(3,3,(T)(sPtr->q*wgp.x[p]*wgp.y[p]),alpha+9*gi,3,vel.x+p,sPtr->Np,dJ,1);
            for(unsigned int dim = 0; dim < 3; ++dim){
                J[3*gi+dim] +=dJ[dim];
            }
        }
    }
}

template<typename T>
void Field2D3V<T>::accumulateM(const std::vector<Species<T,2,3>*> &species, T dt) {
    auto spacings = this->grid.getSpacings();
    T vol = spacings[0]*spacings[1];
    std::fill(Mg,Mg+9*Ng,(T)0.0);
    std::fill(Mgdx,Mgdx+9*Ng,(T)0.0);
    std::fill(Mgdy,Mgdy+9*Ng,(T)0.0);
    std::fill(Mgdxdy,Mgdxdy+9*Ng,(T)0.0);
    std::fill(Mgmdxdy,Mgmdxdy+9*Ng,(T)0.0);
    for(auto sPtr: species){
        const Vector2<unsigned int> g = ((Species2D3V<T>*)sPtr)->getG();
        const Vector2<unsigned int> gp = ((Species2D3V<T>*)sPtr)->getGp();
        const Vector2<T> wg = ((Species2D3V<T>*)sPtr)->getWg();
        const Vector2<T> wgp = ((Species2D3V<T>*)sPtr)->getWgp();
        const T* alpha = ((Species2D3V<T>*)sPtr)->getAlpha();
        T factor = vol * sPtr->q * sPtr->q * dt / (2 * sPtr->m);
        for(unsigned int p = 0; p < sPtr->Np; ++p){
            unsigned int gi = g.x[p]+g.y[p]*Nx;
            unsigned int gdxi = gp.x[p]+g.y[p]*Nx;
            unsigned int gdyi = g.x[p]+gp.y[p]*Nx;
            unsigned int gdxdyi = gp.x[p]+gp.y[p]*Nx;
            T weight;
            //Grid cells onto themselves:
            //0 onto 0, stored in 0's Mg
            weight = wg.x[p]*wg.y[p]*wg.x[p]*wg.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mg[gi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //dx onto dx, stored in dx's Mg
            weight = wgp.x[p]*wg.y[p]*wgp.x[p]*wg.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mg[gdxi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //dy onto dy, stored in dy's Mg
            weight = wg.x[p]*wgp.y[p]*wg.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mg[gdyi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //dxdy onto dxdy, stored in dxdy's Mg
            weight = wgp.x[p]*wgp.y[p]*wgp.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mg[gdxdyi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //Grid cells onto cell+dx:
            //0 onto dx, stored in 0's Mgdx
            weight = wg.x[p]*wg.y[p]*wgp.x[p]*wg.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mgdx[gi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //dy onto dxdy, stored in dy's Mgdx
            weight = wg.x[p]*wgp.y[p]*wgp.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mgdx[gdyi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //Grid cells onto cell+dy:
            //0 onto dy, stored in 0's Mgdy
            weight = wg.x[p]*wg.y[p]*wg.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mgdy[gi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //dx onto dxdy, stored in dx's Mgdy
            weight = wgp.x[p]*wg.y[p]*wgp.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mgdy[gdyi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //Grid cell onto cell + dx + dy:
            //0 onto dx, stored in 0's Mgdxdy
            weight = wg.x[p]*wg.y[p]*wgp.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mgdxdy[gi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
            //Grid cell onto cell - dx + dy:
            //dx onto dy, stored in dx's Mgmdxdy
            weight = wgp.x[p]*wg.y[p]*wg.x[p]*wgp.y[p];
            for(unsigned int comp = 0; comp < 9; ++comp){
                Mgmdxdy[gdxi*9+comp] += alpha[p*9+comp]*weight*factor;
            }
        }
    }
}

template<typename T>
void Field2D3V<T>::solveAndAdvance(T dt) {
    //TODO: system solve

}

template class Field2D3V<float>;
template class Field2D3V<double>;