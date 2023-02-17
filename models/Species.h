#include <fstream>

class Distribution;
class Field;

template <typename T>
class Species{
public:
	const int Nd; //Dimensionality of space
	const int Nv; //Dimensionality of velocity
	
	const int Np;
	const T m;
	const T q;
	
	Species() = delete;
	Species(int Np, T m, T q);
	Species(const Species<T>& other); 
	~Species();

	virtual void initializePositions(const Distribution& dist) = 0;
	virtual void initializeVelocities(const Distribution& dist) = 0;
	virtual void initializeFromFile(const std::ofstream& file) = 0;
	
	virtual void advancePositions(T dt, const Field& field) = 0;
	virtual void advanceVelocities(T dt, const Field& field) = 0;
	virtual void updateAlphaAndWeights(const Field& field) = 0;
	
	virtual void savePositions(std::ofstream& outputFile) const = 0;
	virtual void saveVelocities(std::ofstream& outputFile) const = 0;
		
private:
	virtual void computeAlphas(const Field& field) = 0;
	virtual void computeWeights(const Field& field) = 0;
};