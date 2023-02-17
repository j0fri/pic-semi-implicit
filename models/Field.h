#include <fstream>
#include <vector>

class Species;

template <typename T>
class Field{
public:
	const int Nd; //Dimensionality of space
	const int Nx; //Number of grid points in each direction
	const T Lx; //Field size in each direction
	const T dx; //
	const T c; //Speed of light
	const T e0; //Permittivity
	
	Field() = delete;
	Field(int Nx, T Lx, T c, T e0);
	Field(const Field<T>& other);
	~Field();
	
	virtual void initialize(const std::vector<Species>& species) = 0;
	virtual void advanceField(const std::vector<Species>& species) = 0;
	virtual void save(std::ofstream& outputFile) const = 0;
};