#include <stdexcept>
#include <cmath>


//TODO: ref: Statistical Physics (2nd Edition), F. Mandl, Manchester Physics, John Wiley & Sons, 2008, ISBN 9780471915331
template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::Boltzmann(T m, T Kb, T T0){
    if constexpr (Nd == 0 || Nd > 3){
        throw std::invalid_argument("Nd must be at most 3");
    }
    std::function<T(const std::array<T,Nd>&)> fun;
    if constexpr (Nd == 1){
        T a1 = std::sqrt(m/(2*M_PI*Kb*T0));
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*arr[0]*arr[0]);
        });
    }
    if constexpr (Nd == 2){
        T a1 = 2*M_PI*Kb*T0;
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*(arr[0]*arr[0]+arr[1]*arr[1]));
        });
    }
    if constexpr (Nd == 3){
        T a1 = std::pow(m/(2*M_PI*Kb*T0),1.5);
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*(arr[0]*arr[0]+arr[1]*arr[1]+arr[2]*arr[2]));
        });
    }
    return Distribution<T,Nd>(fun);
}


template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::ShiftedBoltzmann(T m, T Kb, T T0, const std::array<T,Nd>& u0){
    if constexpr (Nd == 0 || Nd > 3){
        throw std::invalid_argument("Nd must be at most 3");
    }
    std::function<T(const std::array<T,Nd>&)> fun;
    if constexpr (Nd == 1){
        T a1 = std::sqrt(m/(2*M_PI*Kb*T0));
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*(arr[0]-u0[0])*(arr[0]-u0[0]));
        });
    }
    if constexpr (Nd == 2){
        T a1 = 2*M_PI*Kb*T0;
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*((arr[0]-u0[0])*(arr[0]-u0[0])+(arr[1]-u0[1])*(arr[1]-u0[1])));
        });
    }
    if constexpr (Nd == 3){
        T a1 = std::pow(m/(2*M_PI*Kb*T0),1.5);
        T a2 = m/(2*Kb*T0);
        fun = std::function([=](const std::array<T,Nd>& arr){
            return a1 * std::exp(-a2*((arr[0]-u0[0])*(arr[0]-u0[0])+(arr[1]-u0[1])*(arr[1]-u0[1])+(arr[2]-u0[2])*(arr[2]-u0[2])));
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


template <typename T, unsigned int Nd>
Distribution<T,Nd> preset_distributions::TopHat(unsigned int Id, T x1, T x2){
    if(Id > Nd){
        throw std::invalid_argument("Id must be smaller than or equal to Nd.");
    }
    auto fun = std::function([=](const std::array<T,Nd>& arr){return arr[Id] >= x1 ? (arr[Id] <= x2 ? (T)1 : (T)0) : (T)0;});
    return Distribution<T,Nd>(fun);
}


template <typename T, unsigned int Nd>
DistributionGrid<T,Nd> preset_distributions::Constant(T val){
    auto fun = std::function([=](const std::array<T,Nd>& arr) {return 1;});
    Grid<T,Nd> grid{};
    for(unsigned int i = 0; i < Nd; ++i){
        grid.dimensions[i] = {val,val,1};
    }
    return DistributionGrid<T,Nd>(fun,grid);

}