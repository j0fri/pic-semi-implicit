#include <fstream>
#include <vector>
#include <array>
#include "Grid.h"

template <typename T, unsigned int Nd, unsigned int Nv>
class Species;

template <typename T, unsigned int Nd, unsigned int Nv>
class Field{
public:
    const Grid<T,Nd> grid; //Grid (including length and spacing)
	const T c; //Speed of light
	const T e0; //Permittivity

	Field() = delete;
	Field(const Grid<T,Nd>& grid, T c, T e0);
	Field(const Field<T,Nd,Nv>& other);
	~Field();
	
	virtual void initialize(const std::vector<Species<T,Nd,Nv>>& species) = 0;
	virtual void advanceField(const std::vector<Species<T,Nd,Nv>>& species) = 0;
	virtual void save(std::ofstream& outputFile) const = 0;
};