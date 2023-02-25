#ifndef PIC_SEMI_IMPLICIT_SPECIES_H
#define PIC_SEMI_IMPLICIT_SPECIES_H

#include <fstream>

#include "Distribution.h"
#include "Config.h"

template <typename T, unsigned int Nd, unsigned int Nv>
class Field;

template <typename T, unsigned int Nd, unsigned int Nv>
class Species{
public:
	const unsigned int Np;
	const T m;
	const T q;
    const Distribution<T,Nd> initialXDist;
    const Distribution<T,Nv> initialVDist;
    const bool initialisePositionFromFile;
    const std::string initialPositionFileName;
    const bool initialiseVelocityFromFile;
    const std::string initialVelocityFileName;
	
	Species() = delete;
	Species(const Config<T,Nd,Nv>::SpeciesConfig& speciesConfig);
	Species(const Species<T,Nd,Nv>& other);
	~Species();

    void initialise();

    virtual void advancePositions(T dt, const Field<T,Nd,Nv>& field) = 0;
    virtual void advanceVelocities(T dt, const Field<T,Nd,Nv>& field) = 0;
    virtual void save(std::ofstream& outputFile, const Config<T,Nd,Nv>::SaveConfig& saveConfig) const = 0;

private:
	virtual void initialisePositions() = 0;
	virtual void initialiseVelocities() = 0;
	virtual void initialisePositions(const std::ifstream& file) = 0;
    virtual void initialiseVelocities(const std::ifstream& file) = 0;

	virtual void updateAlphaAndWeights(const Field<T,Nd,Nv>& field) = 0;

	virtual void computeAlphas(const Field<T,Nd,Nv>& field) = 0;
	virtual void computeWeights(const Field<T,Nd,Nv>& field) = 0;
};

#endif //PIC_SEMI_IMPLICIT_SPECIES_H