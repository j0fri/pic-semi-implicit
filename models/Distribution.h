#ifndef PIC_SEMI_IMPLICIT_DISTRIBUTION_H
#define PIC_SEMI_IMPLICIT_DISTRIBUTION_H

#include <functional>
#include <array>
#include <vector>
#include "Grid.h"


template <typename T, unsigned int Nd>
class Distribution
{
public:
    std::function<T(const std::array<T,Nd>&)> f;
	Distribution() = delete;
	explicit Distribution(const std::function<T(std::array<T,Nd>)>& f);
	std::vector<std::array<T, Nd>> generate(int Np, const Grid<T,Nd>& grid) const;
    Distribution<T,Nd> operator+(const Distribution<T,Nd>& other) const;
    Distribution<T,Nd> operator-(const Distribution<T,Nd>& other) const;
    Distribution<T,Nd> operator*(const Distribution<T,Nd>& other) const;
    Distribution<T,Nd> operator/(const Distribution<T,Nd>& other) const;
    Distribution<T,Nd> operator+=(const Distribution<T,Nd>& other);
    Distribution<T,Nd> operator-=(const Distribution<T,Nd>& other);
    Distribution<T,Nd> operator*=(const Distribution<T,Nd>& other);
    Distribution<T,Nd> operator/=(const Distribution<T,Nd>& other);
private:
	struct Cell {
		std::array<T,Nd> left;
		std::array<T,Nd> right;
		T objectiveNp;
		unsigned int Np;
		std::array<T,Nd> getCentre() const;
	};
	std::vector<Cell> generateMesh(const Grid<T,Nd>& config) const;
    std::vector<std::array<T, Nd>> generate(int Np, const Grid<T,Nd>& grid, T particleIncreaseGenerationFactor) const;

};


#endif // PIC_SEMI_IMPLICIT_DISTRIBUTION_H

