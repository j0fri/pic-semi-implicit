#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <functional>
#include <array>
#include <vector>


template <typename T, unsigned int Nd>
class Distribution
{
public:
    const std::function<T(std::array<T,Nd>)> f;

	struct Config{
		struct Dim{
			T min;
			T max;
			unsigned int Nc;
		};
		std::array<Dim,Nd> dimensions;
        std::array<T,Nd> getSpacings() const;
	};

	Distribution() = delete;
	Distribution(const std::function<T(std::array<T,Nd>)>& f);
	std::vector<std::array<T, Nd>> generate(int Np, const Config& config) const;
private:
	struct Cell {
		std::array<T,Nd> left;
		std::array<T,Nd> right;
		T objectiveNp;
		unsigned int Np;
		std::array<T,Nd> getCentre() const;
	};
	std::vector<Cell> generateMesh(const Config& config) const;
};


#endif // DISTRIBUTION_H

