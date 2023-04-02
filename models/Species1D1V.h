#ifndef PIC_SEMI_IMPLICIT_SPECIES1D1V_H
#define PIC_SEMI_IMPLICIT_SPECIES1D1V_H

#include "Species.h"

template <typename T>
class Species1D1V: public Species<T,1,1> {
private:
    T* x;
    T* v;
    int* g;
    int* gp;
    T* wg;
    T* wgp;
    T* Ep;

public:
    //TODO: add move constructor and assignment operator
    //TODO: actually implement these
    Species1D1V() = delete;
    explicit Species1D1V(const typename Config<T,1,1>::SpeciesConfig& speciesConfig,
                         const typename Config<T,1,1>::BCConfig& bcConfig);
    Species1D1V(const Species1D1V<T>& other);
    ~Species1D1V();

    const T* getV() const;
    const int* getG() const;
    const int* getGp() const;
    const T* getWg() const;
    const T* getWgp() const;
    T getTotalKineticEnergy() const;

    void advancePositions(T dt, const Field<T,1,1>* field);
    void advanceVelocities(T dt, const Field<T,1,1>* field);

    //TODO: implement these
    void savePosition(std::ofstream& outputFile) const;
    void savePositionDistribution(std::ofstream& outputFile) const;
    void saveVelocity(std::ofstream& outputFile) const;
    void saveVelocityDistribution(std::ofstream& outputFile) const;
    void saveEnergy(std::ofstream& outputFile) const;

private:
    void initialisePositions();
    void initialiseVelocities();
    void initialisePositions(const std::ifstream& file);
    void initialiseVelocities(const std::ifstream& file);

    void computeAlphas(const Field<T,1,1>* field, T dt);
    void computeWeights(const Field<T,1,1>* field);
};


#endif //PIC_SEMI_IMPLICIT_SPECIES1D1V_H
