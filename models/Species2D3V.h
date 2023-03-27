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
    Vector2<unsigned int> g; //This gives index of largest cell with cell_coord <= part_coord, where coord is x,y
    Vector2<unsigned int> gp; //This gives index of smallest cell with cell_coord > part_coord, where coord is x,y
    Vector2<unsigned int> gB; //Same but with magnetic field reference
    Vector2<unsigned int> gpB; //Same but with magnetic field reference
    //For non-periodic boundary conditions, particles will have both g and gp in range as in the periodic case,
    //but the weight will be zero if appropriate.
    Vector2<T> wg;
    Vector2<T> wgp;
    Vector2<T> wgB;
    Vector2<T> wgpB;
    Vector3<T> Ep; //Electric field at particle at time n+theta
    Vector3<T> Bp; //Magnetic field at particle at time n
    T* alpha; //Np*9 array, where every alpha is stored in column major format sequentially for every particles

public:
    //TODO: add move constructor and assignment operator
    //TODO: actually implement these
    Species2D3V() = delete;
    explicit Species2D3V(const typename Config<T,2,3>::SpeciesConfig& speciesConfig, const Config<T,2,3>::BCConfig& bcConfig);
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

    void computeAlphas(const Field<T,2,3>* field, T dt);
    void computeWeights(const Field<T,2,3>* field);
    void computeLocalB(const Field<T,2,3>* field);
    void computeLocalE(const Field<T,2,3>* field);
};


#endif //PIC_SEMI_IMPLICIT_SPECIES2D3V_H
