#include "Distribution.h"

#ifndef PIC_SEMI_IMPLICIT_CONFIG_H
#define PIC_SEMI_IMPLICIT_CONFIG_H

#include <istream>

template <typename T, unsigned int Nd, unsigned int Nv>
struct Config{
    struct SpeciesConfig{
        unsigned int Np{};
        T m; //Mass of super-particles
        T q; //Charge of super-particles
        Distribution<T,Nd> xDist;
        Distribution<T,Nv> vDist;
        //TODO: add file for each species
        bool initialisePositionFromFile{};
        std::string initialPositionFileName;
        bool initialiseVelocityFromFile{};
        std::string initialVelocityFileName;
    };
    std::vector<SpeciesConfig> speciesConfig;

    struct FieldConfig{
        Grid<T,Nd> grid;
        T c;
        T e0;
        Distribution<T,Nd> forcedE;
        Distribution<T,Nd> forcedB;
        bool onlyForcedE; //If true, electric field is only forcedE, if false, forcedE is added to PDE sol
        bool onlyForcedB; //If true, magnetic field is only forcedB, if false, forcedB is added to PDE sol
        bool initialiseFromSpecies; //If true, fields are initialised to satisfy divergence laws.
    };
    FieldConfig fieldConfig;

    struct TimeConfig{
        T total;
        T step;
    };
    TimeConfig timeConfig;

    struct SaveConfig{
        bool savePosition{};
        bool savePositionDistribution{};
        bool saveVelocity{};
        bool saveVelocityDistribution{};
        bool saveEnergies{};
        bool saveElectricField{};
        bool saveMagneticField{};
        bool saveVoltage{};
        T saveInterval;
        std::string speciesOutputFileName;
        std::string fieldOutputFileName;
    };
    SaveConfig saveConfig;
};

#endif //PIC_SEMI_IMPLICIT_CONFIG_H
