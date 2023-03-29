#include "Distribution.h"

#ifndef PIC_SEMI_IMPLICIT_CONFIG_H
#define PIC_SEMI_IMPLICIT_CONFIG_H

#include <istream>
#include <optional>

template <typename T, unsigned int Nd, unsigned int Nv>
struct Config{
    struct SpeciesConfig{
        unsigned int Np{};
        T m; //Mass of super-particles
        T q; //Charge of super-particles
        Distribution<T,Nd> xDist;
        Grid<T,Nd> initialXGrid;
        Distribution<T,Nv> vDist;
        Grid<T,Nv> initialVGrid;
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
        std::array<Distribution<T,Nd>,Nv> forcedE;
        std::array<Distribution<T,Nd>,Nv> forcedB;
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
        bool saveSpeciesEnergy{};
        bool saveElectricField{};
        bool saveMagneticField{};
        bool saveFieldEnergy{};
        bool saveVoltage{};
        T saveInterval;
        std::string outputFilesDirectory{"outputs/"};
        std::string outputFilesSubscript{}; //Will be appended to every output file

        //Default file values:
        std::string speciesPositionFileName{"speciesPosition.txt"};
        std::string speciesPositionDistributionFileName{"speciesPositionDistribution.txt"};
        std::string speciesVelocityFileName{"speciesVelocity.txt"};
        std::string speciesVelocityDistributionFileName{"speciesVelocityDistribution.txt"};
        std::string speciesEnergyFileName{"speciesEnergy.txt"};
        std::string electricFieldFileName{"electricField.txt"};
        std::string magneticFieldFileName{"magneticField.txt"};
        std::string fieldEnergyFileName{"fieldEnergy.txt"};
        std::string voltageFileName{"voltage.txt"};
    };
    SaveConfig saveConfig;

    struct BCConfig{
        std::array<bool,Nd> periodic; //for each dimension, if true the boundary will be periodic
        std::array<std::optional<Distribution<T,Nd-1>>, 2*Nd> generators; /*for each dimension, optional generator of
        new particles, order of boundaries is x=minx, x=maxx, y=miny, ..., ignored if periodic in given dimension,
        magnitude/frequency of generator not relevant as created number will correspond to deleted particles*/
    };
    BCConfig bcConfig;

    bool verbose; //Prints time during simulation
};

#endif //PIC_SEMI_IMPLICIT_CONFIG_H
