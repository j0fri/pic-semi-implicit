#include <fstream>

#include "Distribution.h"

template <typename T, unsigned int Nd, unsigned int Nv>
class Field;

template <typename T, unsigned int Nd, unsigned int Nv>
class Species{
public:
	const unsigned int Np;
	const T m;
	const T q;
	
	Species() = delete;
	Species(unsigned int Np, T m, T q);
	Species(const Species<T,Nd,Nv>& other);
	~Species();

	virtual void initializePositions(const Distribution<T,Nd>& dist) = 0;
	virtual void initializeVelocities(const Distribution<T,Nv>& dist) = 0;
	virtual void initializeFromFile(const std::ofstream& file) = 0;
	
	virtual void advancePositions(T dt, const Field<T,Nd,Nv>& field) = 0;
	virtual void advanceVelocities(T dt, const Field<T,Nd,Nv>& field) = 0;
	virtual void updateAlphaAndWeights(const Field<T,Nd,Nv>& field) = 0;
	
	virtual void savePositions(std::ofstream& outputFile) const = 0;
	virtual void saveVelocities(std::ofstream& outputFile) const = 0;
		
private:
	virtual void computeAlphas(const Field<T,Nd,Nv>& field) = 0;
	virtual void computeWeights(const Field<T,Nd,Nv>& field) = 0;
};