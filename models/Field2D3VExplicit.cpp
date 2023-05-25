#include "Field2D3VExplicit.h"


template<typename T>
Field2D3VExplicit<T>::Field2D3VExplicit(const typename Config<T, 2, 3>::FieldConfig &fieldConfig,
                                        const typename Config<T, 2, 3>::BCConfig &bcConfig) :
                                        Field2D3V<T>(fieldConfig, bcConfig) {}

template<typename T>
void Field2D3VExplicit<T>::initialise(const std::vector<Species<T, 2, 3> *> &species, T dt) {
    Field2D3V<T>::initialise(species, dt);

    std::cout << "Solving explicit matrix:" << std::endl;
    solver.analyzePattern(this->A);
    solver.factorize(this->A);
    std::cout << "Explicit matrix solved." << std::endl;
}

template<typename T>
void Field2D3VExplicit<T>::accumulateM(const std::vector<Species<T, 2, 3> *> &species, T dt) {}

template<typename T>
void Field2D3VExplicit<T>::solveAndAdvance(T dt) {
    if(this->onlyForcedE && this->onlyForcedB) {
        //If only forced fields there is no need to solve system
        std::copy(this->field, this->field + 6*this->Ng, this->fieldT);
        return;
    }

    if(this->bcConfig.type==Config<T,2,3>::BC::Periodic){
        this->constructPeriodicC(dt);
    }
    if(this->bcConfig.type==Config<T,2,3>::BC::TwoPlates){
        this->constructTwoPlatesC(dt);
    }

    Eigen::VectorX<T> sol = solver.solve(this->C);

    unsigned int lda = 6*this->Ng;
    unsigned int eq = 0;
    while(eq < lda){
        //TODO: consider if forced terms must be manually added
        if(!this->onlyForcedE){
            for(unsigned int dim = 0; dim < 3; ++dim){
                this->fieldT[eq] = sol[eq];
                this->field[eq] = 2 * sol[eq] - this->field[eq]; //E(n+1) = 2*E(n+1/2)-E(n)
                ++eq;
            }
        }else{
            for(unsigned int dim = 0; dim < 3; ++dim){
                this->fieldT[eq] = this->field[eq];
                ++eq;
            }
        }
        if(!this->onlyForcedB){
            for(unsigned int dim = 0; dim < 3; ++dim){
                this->fieldT[eq] = sol[eq];
                this->field[eq] = 2 * sol[eq] - this->field[eq]; //B(n+1) = 2*B(n+1/2)-B(n)
                ++eq;
            }
        }else{
            for(unsigned int dim = 0; dim < 3; ++dim){
                this->fieldT[eq] = this->field[eq];
                ++eq;
            }
        }
    }
}

template class Field2D3VExplicit<float>;
template class Field2D3VExplicit<double>;