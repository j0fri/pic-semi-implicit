#ifndef PIC_SEMI_IMPLICIT_SPECIES_H
#define PIC_SEMI_IMPLICIT_SPECIES_H

#include <fstream>

#include "Distribution.h"
#include "Config.h"
#include "Field.h"

template <typename T, unsigned int Nd, unsigned int Nv>
class Field;

template <typename T, unsigned int Nd, unsigned int Nv>
class Species{
public:
	const unsigned int Np;
	const T m;
	const T q;
    const Distribution<T,Nd> initialXDist;
    const Grid<T,Nd> initialXGrid;
    const Distribution<T,Nv> initialVDist;
    const Grid<T,Nv> initialVGrid;
    const bool initialisePositionFromFile;
    const std::string initialPositionFileName;
    const bool initialiseVelocityFromFile;
    const std::string initialVelocityFileName;
	
	Species() = delete;
	explicit Species(const Config<T,Nd,Nv>::SpeciesConfig& speciesConfig);
	Species(const Species<T,Nd,Nv>& other);
	virtual ~Species();

    void initialise();

    //TODO: ADD BOUNDARY CONDITIONS
    virtual void advancePositions(T dt, const Field<T,Nd,Nv>* field) = 0;
    virtual void advanceVelocities(T dt, const Field<T,Nd,Nv>* field) = 0;

    virtual void savePosition(std::ofstream& outputFile) const = 0;
    virtual void savePositionDistribution(std::ofstream& outputFile, const Field<T, 1, 1> *field) const = 0;
    virtual void saveVelocity(std::ofstream& outputFile) const = 0;
    virtual void saveVelocityDistribution(std::ofstream& outputFile) const = 0;
    virtual void saveEnergy(std::ofstream& outputFile) const = 0;


private:
	virtual void initialisePositions() = 0;
	virtual void initialiseVelocities() = 0;
	virtual void initialisePositions(const std::ifstream& file) = 0;
    virtual void initialiseVelocities(const std::ifstream& file) = 0;

	virtual void computeAlphas(const Field<T,Nd,Nv>* field) = 0;
	virtual void computeWeights(const Field<T,Nd,Nv>* field) = 0;
};

#endif //PIC_SEMI_IMPLICIT_SPECIES_H