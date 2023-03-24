//
// Created by jf1519 on 22/03/23.
//

#ifndef PIC_SEMI_IMPLICIT_SPECIES2D3V_H
#define PIC_SEMI_IMPLICIT_SPECIES2D3V_H

#include "Species.h"
#include "Config.h"
#include "Vector2.h"
#include "Vector3.h"

template <typename T>
class Species2D3V: public Species<T,2,3> {
private:
    Vector2<T> pos;
    Vector3<T> vel;
    Vector2<unsigned int> g;
    Vector2<unsigned int> gp;
    Vector2<T> wg;
    Vector2<T> wgp;
    Vector3<T> Ep;
    Vector3<T> Bp;

public:
    //TODO: add move constructor and assignment operator
    //TODO: actually implement these
    Species2D3V() = delete;
    explicit Species2D3V(const Config<T,2,3>::SpeciesConfig& speciesConfig);
    Species2D3V(const Species2D3V<T>& other);
    ~Species2D3V();

    const Vector3<T>& getVel() const;
    const Vector2<unsigned int>& getG() const;
    const Vector2<unsigned int>& getGp() const;
    const Vector2<T>& getWg() const;
    const Vector2<T>& getWgp() const;
    T getTotalKineticEnergy() const;

    void advancePositions(T dt, const Field<T,2,3>* field);
    void advanceVelocities(T dt, const Field<T,2,3>* field);

    //TODO: implement these
    void savePosition(std::ofstream& outputFile) const;
    void savePositionDistribution(std::ofstream& outputFile, const Field<T, 2, 3> *field) const;
    void saveVelocity(std::ofstream& outputFile) const;
    void saveVelocityDistribution(std::ofstream& outputFile) const;
    void saveEnergy(std::ofstream& outputFile) const;

private:
    void initialisePositions();
    void initialiseVelocities();
    void initialisePositions(const std::ifstream& file);
    void initialiseVelocities(const std::ifstream& file);

    void computeAlphas(const Field<T,2,3>* field);
    void computeWeights(const Field<T,2,3>* field);
};



#endif //PIC_SEMI_IMPLICIT_SPECIES2D3V_H
