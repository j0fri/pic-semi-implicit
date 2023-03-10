#ifndef PIC_SEMI_IMPLICIT_FIELD_H
#define PIC_SEMI_IMPLICIT_FIELD_H

#include <fstream>
#include <vector>
#include <array>

#include "Grid.h"
#include "Distribution.h"
#include "Config.h"
#include "Species.h"

template <typename T, unsigned int Nd, unsigned int Nv>
class Species;

template <typename T, unsigned int Nd, unsigned int Nv>
class Field{
public:
    const Grid<T,Nd> grid; //Grid (including length and spacing)
	const T c; //Speed of light
	const T e0; //Permittivity
    const Distribution<T,Nd> forcedE;
    const Distribution<T,Nd> forcedB;
    const bool onlyForcedE; //If true, electric field is only forcedE, if false, forcedE is added to PDE sol
    const bool onlyForcedB; //If true, magnetic field is only forcedB, if false, forcedB is added to PDE sol
    const bool initialiseFromSpecies;

	Field() = delete;
	explicit Field(const Config<T,Nd,Nv>::FieldConfig& fieldConfig);
	Field(const Field<T,Nd,Nv>& other);
	virtual ~Field();
	
	virtual void initialise(const std::vector<Species<T,Nd,Nv>*>& species) = 0;
	virtual void advanceField(const std::vector<Species<T,Nd,Nv>*>& species, T dt);

	virtual void saveElectricField(std::ofstream& outputFile) const = 0;
    virtual void saveMagneticField(std::ofstream& outputFile) const = 0;
    virtual void saveEnergy(std::ofstream& outputFile) const = 0;
    virtual void saveVoltage(std::ofstream& outputFile) const = 0;

    //TODO: remove this
    virtual const T* getEt() const = 0;
private:
    virtual void accumulateJ(const std::vector<Species<T,Nd,Nv>*>& species) = 0;
    virtual void accumulateM(const std::vector<Species<T,Nd,Nv>*>& species, T dt) = 0;
    virtual void solveAndAdvance(T dt) = 0;
};


#endif //PIC_SEMI_IMPLICIT_FIELD_H