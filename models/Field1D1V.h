#ifndef PIC_SEMI_IMPLICIT_FIELD1D1V_H
#define PIC_SEMI_IMPLICIT_FIELD1D1V_H

#include "Field.h"

template <typename T>
class Field1D1V: public Field<T,1,1> {
private:
    T* E;
    T* Et;

    T* J;
    T* Mgg;
    T* Mggp;

    T* A; //System matrix. Optimisation: avoid direct matrix formulation and use banded matrix
    T* C;

public:
    unsigned int Nx; //TODO: duplicated, remove, otherwise make const
    T dx; //TODO: duplicated, remove, otherwise make const

    //TODO: add move constructor and assignment operator
    //TODO: actually implement these
    Field1D1V() = delete;
    explicit Field1D1V(const Config<T,1,1>::FieldConfig& fieldConfig);
    Field1D1V(const Field1D1V& other);
    ~Field1D1V();

    void initialise(const std::vector<Species<T,1,1>*>& species);

    //TODO: implement these
    void saveElectricField(std::ofstream& outputFile) const ;
    void saveMagneticField(std::ofstream& outputFile) const;
    void saveEnergy(std::ofstream& outputFile) const;
    void saveVoltage(std::ofstream& outputFile) const;

    const T* getEt() const;
private:
    void accumulateJ(const std::vector<Species<T,1,1>*>& species);
    void accumulateM(const std::vector<Species<T,1,1>*>& species, T dt);
    void solveAndAdvance(T dt);
};


#endif //PIC_SEMI_IMPLICIT_FIELD1D1V_H
