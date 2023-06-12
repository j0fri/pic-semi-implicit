#include "Distribution.h"
#include "DistributionGrid.h"

#ifndef PIC_SEMI_IMPLICIT_CONFIG_H
#define PIC_SEMI_IMPLICIT_CONFIG_H

#include <istream>
#include <optional>

template <typename T, unsigned int Nd, unsigned int Nv>
struct Config{
    struct SpeciesConfig{
        unsigned int Np{};
        T m; //Total mass of species
        T q; //Total charge of species
        Distribution<T,Nd> xDist;
        Grid<T,Nd> initialXGrid;
        Distribution<T,Nv> vDist;
        Grid<T,Nv> initialVGrid; //IMPORTANT: the grid must be appropriate to accurately represent the distribution
        //TODO: add file for each species
        bool initialisePositionFromFile{};
        std::string initialPositionFileName;
        bool initialiseVelocityFromFile{};
        std::string initialVelocityFileName;

        std::optional<DistributionGrid<T,Nd>> bcPositionGenerator; /*optional generator of new particles,
        magnitude/frequency of generator not relevant as created number will correspond to deleted particles*/
        std::optional<DistributionGrid<T,Nv>> bcVelocityGenerator; //optional generator of velocity of new particles
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
        std::optional<T> solverTolerance;
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
        bool saveElectrostaticPotential{};
        bool saveCurrent{};
        bool saveRuntime{true};
        T saveInterval;
        std::string outputFilesDirectory{"outputs/"};
        std::string outputFilesSubscript{".txt"}; //Will be appended to every output file

        //Default file values:
        std::string speciesPositionFileName{"speciesPosition"};
        std::string speciesPositionDistributionFileName{"speciesPositionDistribution"};
        std::string speciesVelocityFileName{"speciesVelocity"};
        std::string speciesVelocityDistributionFileName{"speciesVelocityDistribution"};
        std::string speciesEnergyFileName{"speciesEnergy"};
        std::string electricFieldFileName{"electricField"};
        std::string magneticFieldFileName{"magneticField"};
        std::string fieldEnergyFileName{"fieldEnergy"};
        std::string electrostaticPotentialFileName{"electrostaticPotential"};
        std::string currentFileName{"current"};
        std::string runtimeFileName{"runtime"};

        //If true, saves at all times, ignoring saveInterval
        bool saveAllTimes{false};

        bool saveSolverSteps = false;
        std::string solverStepsFileName{"solverSteps"};
    };
    SaveConfig saveConfig;

    enum BC{Periodic, Diode};
    struct BCConfig{
        BC type;
        struct Diode{
            T voltage;
        };
        std::optional<Diode> diode;
    };
    BCConfig bcConfig;

    bool verbose = false; //Prints time during simulation
    bool outputConfig = false; //Outputs main config values
    bool useExplicitScheme = false; //Use explicit scheme
};

#endif //PIC_SEMI_IMPLICIT_CONFIG_H
