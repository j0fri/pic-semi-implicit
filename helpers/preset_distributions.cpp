#include <stdexcept>
#include <cmath>


//TODO: ref: Statistical Physics (2nd Edition), F. Mandl, Manchester Physics, John Wiley & Sons, 2008, ISBN 9780471915331
template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::Boltzmann(T m, T Kb, T T0){
    if (Nd == 0 || Nd > 3){
        throw std::invalid_argument("Nd must be at most 3");
    }
    std::function<T(const std::array<T,Nd>&)> fun;
    if(Nd == 1){
        T a1 = std::sqrt(m/(2*M_PI*Kb*T0));
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*arr[0]*arr[0]);
        });
    }
    if(Nd == 2){
        T a1 = 2*M_PI*Kb*T0;
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*(arr[0]*arr[0]+arr[1]*arr[1]+arr[2]*arr[2]));
        });
    }
    if(Nd == 3){
        T a1 = std::pow(m/(2*M_PI*Kb*T0),1.5);
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*(arr[0]*arr[0]+arr[1]*arr[1]+arr[2]*arr[2]));
        });
    }
    return Distribution<T,Nd>(fun);
}


template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::Uniform(T a){
    auto fun = std::function([=](const std::array<T,Nd>& arr) {return a;});
    return Distribution<T,Nd>(fun);
}


template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::Step(unsigned int Id, T a, bool reverse){
    if(Id > Nd){
        throw std::invalid_argument("Id must be smaller than or equal to Nd.");
    }
    std::function<T(const std::array<T,Nd>&)> fun;
    if(!reverse){
        fun = std::function([=](const std::array<T,Nd>& arr) {
            return arr[Id]>=a ? 1 : 0;
        });
    }else{
        fun = std::function([=](const std::array<T,Nd>& arr) {
            return arr[Id]<a ? 1 : 0;
        });
    }
    return Distribution<T,Nd>(fun);
}


template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::Sin(unsigned int Id, T a, T b, T c){
    if(Id > Nd){
        throw std::invalid_argument("Id must be smaller than or equal to Nd.");
    }
    auto fun = std::function([=](const std::array<T,Nd>& arr){return a*std::sin(b*arr[Id]+c);});
    return Distribution<T,Nd>(fun);
}