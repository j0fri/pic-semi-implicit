#include "Distribution.h"
#include "../helpers/math_helper.h"

#include <array>
#include <vector>
#include <cmath>
#include <functional>
#include <chrono>
#include <random>
#include <stdexcept>

#define PARTICLE_GENERATION_INCREASE_FACTOR 1.05


template <typename T, unsigned int Nd>
std::array<T,Nd> Distribution<T,Nd>::Cell::getCentre() const{
	std::array<T,Nd> centre{};
	for(unsigned int i = 0; i < Nd; ++i){
		centre[i] = (left[i]+right[i])/2;
	}
	return centre;
}


template <typename T, unsigned int Nd>
Distribution<T,Nd>::Distribution(const std::function<T(std::array<T,Nd>)>& f) : f(f){}


template <typename T, unsigned int Nd>
std::vector<std::array<T, Nd>> Distribution<T,Nd>::generate(int Np, const Grid<T,Nd>& grid) const{
    //TODO: move seeding to Simulation class
    static auto randEngine = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()); //Only seed once

    auto cells = this->generateMesh(grid);
    int Nc = cells.size();

    //Calculate sum of density function to normalize
	T total = 0;
	for (int i = 0; i < Nc; ++i) {
        cells[i].objectiveNp = this->f(cells[i].getCentre());
		if (cells[i].objectiveNp<0){
			throw std::runtime_error("Negative density function at cell centre.");
		}
		total += cells[i].objectiveNp;
	}
	T scaling = 1.0 / total;

    //Allocate particles to each cell
	int totalParticles = 0;
	for (int i = 0; i < Nc; ++i) {
        cells[i].objectiveNp *= Np * scaling * PARTICLE_GENERATION_INCREASE_FACTOR;
        cells[i].Np = std::floor(cells[i].objectiveNp); //Always add at least floor(Npdes)
        T remainder = cells[i].objectiveNp - cells[i].Np;
        if((double)std::rand()/RAND_MAX < remainder){
            ++cells[i].Np; //Add 1 particle with remainder (decimal of Npdes) probability
        }
		totalParticles += cells[i].Np;
	}

    //Check enough particles
    if(totalParticles < Np){
        throw std::runtime_error("Generated particles less than Np. Consider increasing "
                                 "PARTICLE_GENERATION_INCREASE_FACTOR.");
    }

    std::vector<std::array<T,Nd>> out(totalParticles);
    int particleCounter = 0;
	//Distribute particles inside each cell
	for(int i = 0; i < Nc; ++i) {
        const auto& left = cells[i].left;
        const auto& right = cells[i].right;
        for(unsigned int j = 0; j < cells[i].Np; ++j){
            out[particleCounter] = std::array<T,Nd>();
            for(unsigned int k = 0; k < Nd; ++k){
                out[particleCounter][k] = (double)std::rand()/RAND_MAX*(right[k]-left[k])+left[k];
            }
            ++particleCounter;
        }
	}

    //Randomly remove excess particles
    std::shuffle(out.begin(),out.end(),randEngine);
    out.resize(Np);

    return out;
}


template <typename T, unsigned int Nd>
std::vector<typename Distribution<T,Nd>::Cell> Distribution<T,Nd>::generateMesh(const Grid<T,Nd>& grid) const {
    std::array<T,Nd> spacings = grid.getSpacings();
    std::vector<std::vector<T>> linspaces{Nd};
    for(unsigned int i = 0; i < Nd; ++i){
        const auto& dim = grid.dimensions[i];
        linspaces[i] = math_helper::linspace(dim.min, dim.max-spacings[i], dim.Nc);
    }
    std::vector<std::vector<T>> gridPoints = math_helper::repeatedCartesian(linspaces);

    std::vector<Cell> cells{gridPoints.size()};
    for(int i = 0; i < (int)gridPoints.size(); ++i){
        std::array<T,Nd> left{};
        std::array<T,Nd> right{};
        for(unsigned int j = 0; j < Nd; ++j){
            left[j] = gridPoints[i][j];
            right[j] = gridPoints[i][j] + spacings[j];
        }
        cells[i] = {left, right, 0, 0};
    }
    return cells;
}


template<typename T, unsigned int Nd>
Distribution<T, Nd> Distribution<T, Nd>::operator+(const Distribution<T, Nd> &other) const {
    std::function<T(std::array<T,Nd>)> newFun([=,*this](const std::array<T,Nd>& arr){return this->f(arr) + other.f(arr);});
    return Distribution<T, Nd>(newFun);
}


template<typename T, unsigned int Nd>
Distribution<T, Nd> Distribution<T, Nd>::operator-(const Distribution<T, Nd> &other) const {
    std::function<T(std::array<T,Nd>)> newFun([=,*this](const std::array<T,Nd>& arr){return this->f(arr) - other.f(arr);});
    return Distribution<T, Nd>(newFun);
}


template<typename T, unsigned int Nd>
Distribution<T, Nd> Distribution<T, Nd>::operator*(const Distribution<T, Nd> &other) const {
    std::function<T(std::array<T,Nd>)> newFun([=,*this](const std::array<T,Nd>& arr){return this->f(arr) * other.f(arr);});
    return Distribution<T, Nd>(newFun);
}


template<typename T, unsigned int Nd>
Distribution<T, Nd> Distribution<T, Nd>::operator/(const Distribution<T, Nd> &other) const {
    std::function<T(std::array<T,Nd>)> newFun([=,*this](const std::array<T,Nd>& arr){return this->f(arr) / other.f(arr);});
    return Distribution<T, Nd>(newFun);
}


template class Distribution<double,1>;
template class Distribution<double,2>;
template class Distribution<double,3>;
template class Distribution<float,1>;
template class Distribution<float,2>;
template class Distribution<float,3>;