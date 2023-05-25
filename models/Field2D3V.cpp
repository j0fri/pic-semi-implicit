#include "Field2D3V.h"
#include "Species2D3V.h"
#include "../helpers/math_helper.h"
#include "../helpers/output_helper.h"


template<typename T>
Field2D3V<T>::Field2D3V(const typename Config<T,2,3>::FieldConfig &fieldConfig,
                        const typename Config<T,2,3>::BCConfig& bcConfig):
                        Field<T,2,3>(fieldConfig, bcConfig),
                        Ng{fieldConfig.grid.dimensions[0].Nc*fieldConfig.grid.dimensions[1].Nc},
                        Nx{fieldConfig.grid.dimensions[0].Nc}, Ny{fieldConfig.grid.dimensions[1].Nc},
                        A{SpMat(6*Ng, 6*Ng)}, Am{SpMat(6*Ng, 6*Ng)}, Ac{SpMat(6*Ng, 6*Ng)},
                        C{Eigen::VectorX<T>(6*Ng)}, Esolver{}, Ec{Ng} {
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
}


template<typename T>
void Field2D3V<T>::initialise(const std::vector<Species<T, 2, 3> *> &species, T dt) {
    if(this->initialiseFromSpecies){
        std::cout << "WARNING: when initialising field from species, if periodic boundary conditions are present, ensure the "
                     "charge distribution results in a periodic field." << std::endl;
    }
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
    }

    if(this->bcConfig.type==Config<T,2,3>::BC::Periodic){
        this->initialisePeriodicElectrostaticPotentialSystem();
        this->initialisePeriodicField(species, dt);
        this->initialisePeriodicA(dt);
        return;
    }
    if(this->bcConfig.type==Config<T,2,3>::BC::TwoPlates){
        this->initialiseTwoPlateElectrostaticPotentialSystem();
        this->initialiseTwoPlatesField(species, dt);
        this->initialiseTwoPlatesA(dt);
        return;
    }

    throw std::invalid_argument("Boundaries in a 2D3V field must be either all periodic or all perfect conductors.");
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
    T e0 = this->e0;
    T mu0 = (T)1/(this->c*this->c*this->e0);
    T totalE = 0;
    T totalB = 0;
    for(unsigned int gi = 0; gi < this->Ng; ++gi){
        //Electric field terms
        totalE += field[6*gi]*field[6*gi];
        totalE += field[6*gi+1]*field[6*gi+1];
        totalE += field[6*gi+2]*field[6*gi+2];
        //Magnetic field terms
        totalB += field[6*gi+3]*field[6*gi+3];
        totalB += field[6*gi+4]*field[6*gi+4];
        totalB += field[6*gi+5]*field[6*gi+5];
    }
    T dv = this->grid.getSpacings()[0] * this->grid.getSpacings()[1];
    T totalEnergy = totalE*dv*e0/(8*M_PI) + totalB*dv*mu0/(8*M_PI); //Electromagnetic field energy in Gauss units
    outputFile << totalEnergy << std::endl;
}


template<typename T>
void Field2D3V<T>::saveElectrostaticPotential(std::ofstream &outputFile, const std::vector<Species<T,2,3>*> &species) const {
    std::unique_ptr<const T> phi = this->getElectrostaticPotential(species);
    output_helper::outputColMajorMatrix(&*phi,Nx,Ny,Nx,1,outputFile);
}


template<typename T>
void Field2D3V<T>::saveCurrent(std::ofstream &outputFile) const {
    output_helper::outputColMajorMatrix(J,Nx,Ny,Nx*3,3,outputFile);
    output_helper::outputColMajorMatrix(J+1,Nx,Ny,Nx*3,3,outputFile);
    output_helper::outputColMajorMatrix(J+2,Nx,Ny,Nx*3,3,outputFile);
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
    if(this->bcConfig.type==Config<T,2,3>::BC::Periodic){
        return this->getPeriodicElectrostaticPotential(species);
    }
    if(this->bcConfig.type==Config<T,2,3>::BC::TwoPlates){
        return this->getPeriodicElectrostaticPotential(species);
    }

    throw std::invalid_argument("Unsupported boundary conditions in Field2D3V getElectrostaticPotential.");
}


template<typename T>
void Field2D3V<T>::accumulateJ(const std::vector<Species<T,2,3>*> &species) {
    std::fill(J,J+3*Ng,(T)0.0);
    T dV = this->grid.getSpacings()[0] * this->grid.getSpacings()[1];
    for(auto sPtr: species){
        T dJ[3];
        T pVel[3];
        const Vector2<unsigned int> g = ((Species2D3V<T>*)sPtr)->getG();
        const Vector2<unsigned int> gp = ((Species2D3V<T>*)sPtr)->getGp();
        const Vector2<T> wg = ((Species2D3V<T>*)sPtr)->getWg();
        const Vector2<T> wgp = ((Species2D3V<T>*)sPtr)->getWgp();
        const Vector3<T> vel = ((Species2D3V<T>*)sPtr)->getVel();
        const T* alpha = ((Species2D3V<T>*)sPtr)->getAlpha();
        const T* alphaPtr;
        T* jPtr;
        T* dJPtr;
        for(unsigned int p = 0; p < sPtr->Np; ++p){
            pVel[0] = vel.x[p];
            pVel[1] = vel.y[p];
            pVel[2] = vel.z[p];
            //Bottom-left cell
            unsigned int gi = g.x[p]+g.y[p]*Nx;
            alphaPtr = alpha+9*gi;
            math_helper::gemv(3,3,(T)(sPtr->q*wg.x[p]*wg.y[p]),alphaPtr,3,pVel,1,dJ,1);
            jPtr = J + 3*gi-1;
            dJPtr = dJ-1;
            for(unsigned int dim = 0; dim < 3; ++dim){
                *(++jPtr) += *(++dJPtr);
            }
            //Bottom-right cell (dx)
            gi = gp.x[p]+g.y[p]*Nx;
            alphaPtr = alpha+9*gi;
            math_helper::gemv(3,3,(T)(sPtr->q*wgp.x[p]*wg.y[p]),alphaPtr,3,pVel,1,dJ,1);
            jPtr = J + 3*gi-1;
            dJPtr = dJ-1;
            for(unsigned int dim = 0; dim < 3; ++dim){
                *(++jPtr) += *(++dJPtr);
            }
            //Top-left cell (dy)
            gi = g.x[p]+gp.y[p]*Nx;
            alphaPtr = alpha+9*gi;
            math_helper::gemv(3,3,(T)(sPtr->q*wg.x[p]*wgp.y[p]),alphaPtr,3,pVel,1,dJ,1);
            jPtr = J + 3*gi-1;
            dJPtr = dJ-1;
            for(unsigned int dim = 0; dim < 3; ++dim){
                *(++jPtr) += *(++dJPtr);
            }
            //Top-right cell (dx, dy)
            gi = gp.x[p]+gp.y[p]*Nx;
            alphaPtr = alpha+9*gi;
            math_helper::gemv(3,3,(T)(sPtr->q*wgp.x[p]*wgp.y[p]),alphaPtr,3,pVel,1,dJ,1);
            jPtr = J + 3*gi-1;
            dJPtr = dJ-1;
            for(unsigned int dim = 0; dim < 3; ++dim){
                *(++jPtr) += *(++dJPtr);
            }
        }
    }
    T idV = (T)1/dV;
    for(unsigned int i = 0; i < 3*Ng; ++i){
        J[i]*=idV;
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
        T factor = sPtr->q * sPtr->q * dt / (2 * sPtr->m * vol);
        for(unsigned int p = 0; p < sPtr->Np; ++p){
            unsigned int gi9 = 9*(g.x[p]+g.y[p]*Nx);
            unsigned int gdxi9 = 9*(gp.x[p]+g.y[p]*Nx);
            unsigned int gdyi9 = 9*(g.x[p]+gp.y[p]*Nx);
            unsigned int gdxdyi9 = 9*(gp.x[p]+gp.y[p]*Nx);
            T weight;
            T weightFactor;
            T* Mptr;
            const T* alphaPtr;
            const T* const alphaPtrConst = alpha+p*9-1;
            //Grid cells onto themselves:
            //0 onto 0, stored in 0's Mg
            weight = wg.x[p]*wg.y[p]*wg.x[p]*wg.y[p];
            weightFactor = weight*factor;
            Mptr = Mg + gi9 - 1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //dx onto dx, stored in dx's Mg
            weight = wgp.x[p]*wg.y[p]*wgp.x[p]*wg.y[p];
            weightFactor = weight*factor;
            Mptr = Mg+gdxi9-1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //dy onto dy, stored in dy's Mg
            weight = wg.x[p]*wgp.y[p]*wg.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mg+gdyi9-1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //dxdy onto dxdy, stored in dxdy's Mg
            weight = wgp.x[p]*wgp.y[p]*wgp.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mg+gdxdyi9-1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //Grid cells onto cell+dx:
            //0 onto dx, stored in 0's Mgdx
            weight = wg.x[p]*wg.y[p]*wgp.x[p]*wg.y[p];
            weightFactor = weight*factor;
            Mptr = Mgdx + gi9 - 1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //dy onto dxdy, stored in dy's Mgdx
            weight = wg.x[p]*wgp.y[p]*wgp.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mgdx+gdyi9-1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //Grid cells onto cell+dy:
            //0 onto dy, stored in 0's Mgdy
            weight = wg.x[p]*wg.y[p]*wg.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mgdy + gi9 - 1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //dx onto dxdy, stored in dx's Mgdy
            weight = wgp.x[p]*wg.y[p]*wgp.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mgdy+gdyi9-1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //Grid cell onto cell + dx + dy:
            //0 onto dx, stored in 0's Mgdxdy
            weight = wg.x[p]*wg.y[p]*wgp.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mgdxdy + gi9 - 1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
            //Grid cell onto cell - dx + dy:
            //dx onto dy, stored in dx's Mgmdxdy
            weight = wgp.x[p]*wg.y[p]*wg.x[p]*wgp.y[p];
            weightFactor = weight*factor;
            Mptr = Mgmdxdy+gdxi9-1;
            alphaPtr = alphaPtrConst;
            for(unsigned int comp = 0; comp < 9; ++comp){
                *(++Mptr) += *(++alphaPtr)*weightFactor;
            }
        }
    }
}


template<typename T>
void Field2D3V<T>::solveAndAdvance(T dt) {
    if(this->onlyForcedE && this->onlyForcedB) {
        //If only forced fields there is no need to solve system
        std::copy(field, field + 6*Ng, fieldT);
        return;
    }

    if(this->bcConfig.type==Config<T,2,3>::BC::Periodic){
        this->constructPeriodicAc(dt);
        this->constructPeriodicC(dt);
    }
    if(this->bcConfig.type==Config<T,2,3>::BC::TwoPlates){
        this->constructTwoPlatesAc(dt);
        this->constructTwoPlatesC(dt);
    }

    Eigen::LeastSquaresConjugateGradient<SpMat> lscg;
    lscg.setTolerance((T)1e-8);
    lscg.compute(Ac);
    Eigen::VectorX<T> sol = lscg.solve(C);

    unsigned int lda = 6*Ng;
    unsigned int eq = 0;
    while(eq < lda){
        //TODO: consider if forced terms must be manually added
        if(!this->onlyForcedE){
            for(unsigned int dim = 0; dim < 3; ++dim){
                fieldT[eq] = sol[eq];
                field[eq] = 2 * sol[eq] - field[eq]; //E(n+1) = 2*E(n+1/2)-E(n)
                ++eq;
            }
        }else{
            for(unsigned int dim = 0; dim < 3; ++dim){
                fieldT[eq] = field[eq];
                ++eq;
            }
        }
        if(!this->onlyForcedB){
            for(unsigned int dim = 0; dim < 3; ++dim){
                fieldT[eq] = sol[eq];
                field[eq] = 2 * sol[eq] - field[eq]; //B(n+1) = 2*B(n+1/2)-B(n)
                ++eq;
            }
        }else{
            for(unsigned int dim = 0; dim < 3; ++dim){
                fieldT[eq] = field[eq];
                ++eq;
            }
        }
    }
}


template<typename T>
void Field2D3V<T>::initialisePeriodicA(T dt) {
    unsigned int eq = 0;
    unsigned int xi = 0;
    unsigned int yi = 0;
    unsigned int gi, gdxi, gdxdyi, gdyi, gmdxi, gmdxmdyi, gmdyi;
    T dx = this->grid.getSpacings()[0];
    T dy = this->grid.getSpacings()[1];
    T cx = (this->c*dt)/(4*dx);
    T cy = (this->c*dt)/(4*dy);

    //Construct system matrix:
    typedef Eigen::Triplet<T> Tri;
    std::vector<Tri> tripletList;

    while(eq < 6*Ng){
        //Indices:
        gi = xi + Nx*yi;
        gdxi = (xi+1)%Nx + Nx*yi;
        gdxdyi = (xi+1)%Nx + Nx*((yi+1)%Ny);
        gdyi = xi + Nx*((yi+1)%Ny);
        gmdxi = (xi-1+Nx)%Nx + Nx*yi;
        gmdxmdyi = (xi-1+Nx)%Nx + Nx*((yi-1+Ny)%Ny);
        gmdyi = xi + Nx*((yi-1+Ny)%Ny);

        //Electric field equations:
        //x-equation:
        tripletList.emplace_back(eq, 6*gi, (T)1);
        tripletList.emplace_back(eq, 6*gmdxi+5, -cy);
        tripletList.emplace_back(eq, 6*gi+5, -cy);
        tripletList.emplace_back(eq, 6*gmdxmdyi+5, cy);
        tripletList.emplace_back(eq, 6*gmdyi+5, cy);
        ++eq;

        //y-equation:
        tripletList.emplace_back(eq, 6*gi+1, (T)1);
        tripletList.emplace_back(eq, 6*gmdyi+5, cx);
        tripletList.emplace_back(eq, 6*gi+5, cx);
        tripletList.emplace_back(eq, 6*gmdxmdyi+5, -cx);
        tripletList.emplace_back(eq, 6*gmdxi+5, -cx);
        ++eq;

        //z-equation:
        tripletList.emplace_back(eq, 6*gi+2, (T)1);
        tripletList.emplace_back(eq, 6*gmdyi+4, -cx);
        tripletList.emplace_back(eq, 6*gi+4, -cx);
        tripletList.emplace_back(eq, 6*gmdxmdyi+4, cx);
        tripletList.emplace_back(eq, 6*gmdxi+4, cx);
        tripletList.emplace_back(eq, 6*gmdxi+3, cy);
        tripletList.emplace_back(eq, 6*gi+3, cy);
        tripletList.emplace_back(eq, 6*gmdxmdyi+3, -cy);
        tripletList.emplace_back(eq, 6*gmdyi+3, -cy);
        ++eq;

        //Magnetic field equations:
        //x-equation:
        tripletList.emplace_back(eq, 6*gi+3, (T)1);
        tripletList.emplace_back(eq, 6*gdyi+2, cy);
        tripletList.emplace_back(eq, 6*gdxdyi+2, cy);
        tripletList.emplace_back(eq, 6*gi+2, -cy);
        tripletList.emplace_back(eq, 6*gdxi+2, -cy);
        ++eq;

        //y-equation:
        tripletList.emplace_back(eq, 6*gi+4, (T)1);
        tripletList.emplace_back(eq, 6*gdxi+2, -cx);
        tripletList.emplace_back(eq, 6*gdxdyi+2, -cx);
        tripletList.emplace_back(eq, 6*gi+2, cx);
        tripletList.emplace_back(eq, 6*gdyi+2, cx);
        ++eq;

        //z-equation:
        tripletList.emplace_back(eq, 6*gi+5, (T)1);
        tripletList.emplace_back(eq, 6*gdxi+1, cx);
        tripletList.emplace_back(eq, 6*gdxdyi+1, cx);
        tripletList.emplace_back(eq, 6*gi+1, -cx);
        tripletList.emplace_back(eq, 6*gdyi+1, -cx);
        tripletList.emplace_back(eq, 6*gdyi, -cy);
        tripletList.emplace_back(eq, 6*gdxdyi, -cy);
        tripletList.emplace_back(eq, 6*gi, cy);
        tripletList.emplace_back(eq, 6*gdxi, cy);
        ++eq;

        //Advance indices:
        ++xi;
        if(xi == Nx){
            xi = 0;
            ++yi;
        }
    }

    A = SpMat(6*Ng, 6*Ng);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}


template<typename T>
void Field2D3V<T>::initialiseTwoPlatesA(T dt) {
    //TODO
}


template<typename T>
void Field2D3V<T>::initialisePeriodicField(const std::vector<Species<T,2,3>*>& species, T dt) {
    T dx = this->grid.getSpacings()[0];
    T dy = this->grid.getSpacings()[1];

    std::fill(field, field+6*Ng,(T)0.0);

    std::unique_ptr<const T> phi = this->getElectrostaticPotential(species);
    T minx = this->grid.dimensions[0].min;
    T miny = this->grid.dimensions[1].min;

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
}


template<typename T>
void Field2D3V<T>::initialiseTwoPlatesField(const std::vector<Species<T,2,3>*>& species, T dt) {
    //TODO
}


template<typename T>
void Field2D3V<T>::constructPeriodicAc(T dt) {
    //Construct system matrix:
    typedef Eigen::Triplet<T> Tri;
    std::vector<Tri> tripletList;

    unsigned int lda = 6*Ng;
    unsigned int gi, gdxi, gdxdyi, gdyi, gmdxdyi, gmdxi, gmdxmdyi, gmdyi, gdxmdyi;
    T c1 = 2*M_PI*dt;
    unsigned int eq = 0;
    unsigned int xi = 0;
    unsigned int yi = 0;
    while(eq < lda) {
        //Indices:
        gi = xi + Nx*yi;
        gdxi = (xi+1)%Nx + Nx*yi;
        gdxdyi = (xi+1)%Nx + Nx*((yi+1)%Ny);
        gdyi = xi + Nx*((yi+1)%Ny);
        gmdxdyi = (xi-1+Nx)%Nx + Nx*((yi+1)%Ny);
        gmdxi = (xi-1+Nx)%Nx + Nx*yi;
        gmdxmdyi = (xi-1+Nx)%Nx + Nx*((yi-1+Ny)%Ny);
        gmdyi = xi + Nx*((yi-1+Ny)%Ny);
        gdxmdyi = (xi+1)%Nx + Nx*((yi-1+Ny)%Ny);

        //gi term: add Mg term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gi+col, c1 * Mg[9*gi + 3*col + row]);
            }

        }
        //gdxi term: add Mgdx term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdxi+col, c1 * Mgdx[9*gi + 3*col + row]);
            }
        }
        //gdxdyi term: add Mgdxdy term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdxdyi+col, c1 * Mgdxdy[9*gi + 3*col + row]);
            }
        }
        //gdyi term: add Mgdy term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdyi+col, c1 * Mgdy[9*gi + 3*col + row]);
            }
        }
        //gmdxdyi term: add Mgmdxdy
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdxdyi+col, c1 * Mgmdxdy[9*gi + 3*col + row]);
            }
        }
        //gmdxi term: add Mgdx term at gmdxi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdxi+col, c1 * Mgdx[9*gmdxi + 3*col + row]);
            }
        }
        //gmdxmdyi term: add Mgdxdy term at gmdxmdyi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdxmdyi+col, c1 * Mgdxdy[9*gmdxmdyi + 3*col + row]);
            }
        }
        //gmdyi term: add Mgdy term at gmdyi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdyi+col, c1 * Mgdy[9*gmdyi + 3*col + row]);
            }
        }
        //gdxmdyi term: add Mgmdxdy term at gdxmdyi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdxmdyi+col, c1 * Mgmdxdy[9*gdxmdyi + 3*col + row]);
            }
        }
        eq += 6;
        //Advance indices:
        ++xi;
        if(xi == Nx){
            xi = 0;
            ++yi;
        }
    }

    Am = SpMat(6*Ng,6*Ng);
    Am.setFromTriplets(tripletList.begin(), tripletList.end());
    Ac = A;
    Ac += Am;
}


template<typename T>
void Field2D3V<T>::constructTwoPlatesAc(T dt) {
    typedef Eigen::Triplet<T> Tri;
    std::vector<Tri> tripletList;

    unsigned int lda = 6*Ng;
    unsigned int gi, gdxi, gdxdyi, gdyi, gmdxdyi, gmdxi, gmdxmdyi, gmdyi, gdxmdyi;
    T c1 = 2*M_PI*dt;
    unsigned int eq = 0;
    unsigned int xi = 0;
    unsigned int yi = 0;
    while(eq < lda) {
        //Indices:
        gi = xi + Nx*yi;
        gdxi = (xi+1)%Nx + Nx*yi;
        gdxdyi = (xi+1)%Nx + Nx*((yi+1)%Ny);
        gdyi = xi + Nx*((yi+1)%Ny);
        gmdxdyi = (xi-1+Nx)%Nx + Nx*((yi+1)%Ny);
        gmdxi = (xi-1+Nx)%Nx + Nx*yi;
        gmdxmdyi = (xi-1+Nx)%Nx + Nx*((yi-1+Ny)%Ny);
        gmdyi = xi + Nx*((yi-1+Ny)%Ny);
        gdxmdyi = (xi+1)%Nx + Nx*((yi-1+Ny)%Ny);

        //gi term: add Mg term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gi+col, c1 * Mg[9*gi + 3*col + row]);
            }

        }
        //gdxi term: add Mgdx term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdxi+col, c1 * Mgdx[9*gi + 3*col + row]);
            }
        }
        //gdxdyi term: add Mgdxdy term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdxdyi+col, c1 * Mgdxdy[9*gi + 3*col + row]);
            }
        }
        //gdyi term: add Mgdy term at gi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdyi+col, c1 * Mgdy[9*gi + 3*col + row]);
            }
        }
        //gmdxdyi term: add Mgmdxdy
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdxdyi+col, c1 * Mgmdxdy[9*gi + 3*col + row]);
            }
        }
        //gmdxi term: add Mgdx term at gmdxi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdxi+col, c1 * Mgdx[9*gmdxi + 3*col + row]);
            }
        }
        //gmdxmdyi term: add Mgdxdy term at gmdxmdyi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdxmdyi+col, c1 * Mgdxdy[9*gmdxmdyi + 3*col + row]);
            }
        }
        //gmdyi term: add Mgdy term at gmdyi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gmdyi+col, c1 * Mgdy[9*gmdyi + 3*col + row]);
            }
        }
        //gdxmdyi term: add Mgmdxdy term at gdxmdyi
        for(unsigned int col = 0; col < 3; ++col){
            for(unsigned int row = 0; row < 3; ++row){
                tripletList.emplace_back(eq+row,6*gdxmdyi+col, c1 * Mgmdxdy[9*gdxmdyi + 3*col + row]);
            }
        }
        eq += 6;
        //Advance indices:
        ++xi;
        if(xi == Nx){
            xi = 0;
            ++yi;
        }
    }

    Am = SpMat(6*Ng,6*Ng);
    Am.setFromTriplets(tripletList.begin(), tripletList.end());
    Ac = A;
    Ac += Am;
}


template<typename T>
void Field2D3V<T>::constructPeriodicC(T dt) {
    unsigned int lda = 6*Ng;
    T c1 = 2*M_PI*dt;
    unsigned int eq, xi, yi, gi;
    eq = 0;
    xi = 0;
    yi = 0;
    while(eq < lda) {
        //Indices:
        gi = xi + Nx*yi;

        //Electric field equations:
        C[eq++] = field[6 * gi] - c1 * J[3 * gi];
        C[eq++] = field[6 * gi + 1] - c1 * J[3 * gi + 1];
        C[eq++] = field[6 * gi + 2] - c1 * J[3 * gi + 2];

        //Magnetic field equations:
        C[eq++] = field[6 * gi + 3];
        C[eq++] = field[6 * gi + 4];
        C[eq++] = field[6 * gi + 5];

        //Advance indices:
        ++xi;
        if(xi == Nx){
            xi = 0;
            ++yi;
        }
    }
}


template<typename T>
void Field2D3V<T>::constructTwoPlatesC(T dt) {
    //TODO
}


template<typename T>
std::unique_ptr<const T>
Field2D3V<T>::getPeriodicElectrostaticPotential(const std::vector<Species<T, 2, 3> *> &species) const {
    T dx = this->grid.getSpacings()[0];
    T dy = this->grid.getSpacings()[1];

    //Construct system vector: -4 pi rho at every E cell.
    std::fill(Ec.begin(), Ec.end(), (T)0);
    for(auto sPtr: species){
        const Vector2<unsigned int> g = ((Species2D3V<T>*)sPtr)->getG();
        const Vector2<unsigned int> gp = ((Species2D3V<T>*)sPtr)->getGp();
        const Vector2<T> Wg = ((Species2D3V<T>*)sPtr)->getWg();
        const Vector2<T> Wgp = ((Species2D3V<T>*)sPtr)->getWgp();
        T qScaled = -4*M_PI*sPtr->q/(dx*dy); //Gauss units
        for(unsigned int p = 0; p < sPtr->Np; ++p){
            Ec[g.x[p]+this->Nx*g.y[p]] += qScaled*Wg.x[p]*Wg.y[p];
            Ec[g.x[p]+this->Nx*gp.y[p]] += qScaled*Wg.x[p]*Wgp.y[p];
            Ec[gp.x[p]+this->Nx*g.y[p]] += qScaled*Wgp.x[p]*Wg.y[p];
            Ec[gp.x[p]+this->Nx*gp.y[p]] += qScaled*Wgp.x[p]*Wgp.y[p];
        }
    }
    Ec[0]=0; //For potential at 0=0;

    Eigen::VectorX<T> sol = Esolver.solve(Ec);

    T* sol_out = new T[Ng];
    std::copy(sol.begin(), sol.end(), sol_out);
    return std::unique_ptr<const T>(sol_out);
}

template<typename T>
std::unique_ptr<const T>
Field2D3V<T>::getTwoPlatesElectrostaticPotential(const std::vector<Species<T, 2, 3> *> &species) const {
    throw std::invalid_argument("Electrostatic potential not implemented for two-plates.");
}

template<typename T>
void Field2D3V<T>::initialisePeriodicElectrostaticPotentialSystem() {
    T dx = this->grid.getSpacings()[0];
    T dy = this->grid.getSpacings()[1];

    typedef Eigen::Triplet<T> Tri;
    std::vector<Tri> tripletList;

    //Populate system matrix:
    unsigned int eq;
    unsigned int gi, gdxi, gmdxi, gdyi, gmdyi;
    T c1 = (T)1/(dx*dx);
    T c2 = (T)1/(dy*dy);
    for(unsigned int j = 0; j < Ny; ++j){
        for(unsigned int i = 0; i < Nx; ++i){
            eq = j*Nx+i;
            if(eq == 0){
                continue; //This will be the fixed point
            }
            gi = j*Nx+i;
            gdxi = j*Nx+((i+1)%Nx);
            gmdxi = j*Nx+((i-1+Nx)%Nx);
            gdyi = ((j+1)%Ny)*Nx+i;
            gmdyi = ((j-1+Ny)%Ny)*Nx+i;

            //Poisson equation:
            tripletList.emplace_back(eq,gdxi,c1);
            tripletList.emplace_back(eq,gmdxi,c1);
            tripletList.emplace_back(eq,gdyi,c2);
            tripletList.emplace_back(eq,gmdyi,c2);
            tripletList.emplace_back(eq,gi,-2*c1-2*c2);
        }
    }

    //Potential is 0 at 0,0 to avoid matrix rank-deficiency
    tripletList.emplace_back(0,0,1);

    SpMat EA = SpMat(Ng,Ng);
    EA.setFromTriplets(tripletList.begin(), tripletList.end());
    EA.makeCompressed();

    Esolver.analyzePattern(EA);
    Esolver.factorize(EA);
}

template<typename T>
void Field2D3V<T>::initialiseTwoPlateElectrostaticPotentialSystem() {
    throw std::invalid_argument("Electrostatic potential not implemented for two-plates.");
}


template class Field2D3V<float>;
template class Field2D3V<double>;