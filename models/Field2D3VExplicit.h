#ifndef PIC_SEMI_IMPLICIT_FIELD2D3VEXPLICIT_H
#define PIC_SEMI_IMPLICIT_FIELD2D3VEXPLICIT_H

#include "Field2D3V.h"

template<typename T>
class Field2D3VExplicit : public Field2D3V<T>{
    typedef Eigen::SparseMatrix<T> SpMat;
protected:
    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver; //System solver
public:
    Field2D3VExplicit() = delete;
    explicit Field2D3VExplicit(const typename Config<T,2,3>::FieldConfig& fieldConfig,
                       const typename Config<T,2,3>::BCConfig& bcConfig);
    Field2D3VExplicit(const Field2D3VExplicit& other) = delete;
    ~Field2D3VExplicit() = default;

    void initialise(const std::vector<Species<T,2,3>*>& species, T dt);
protected:
    void accumulateM(const std::vector<Species<T,2,3>*>& species, T dt);
    void solveAndAdvance(T dt);
};


#endif //PIC_SEMI_IMPLICIT_FIELD2D3VEXPLICIT_H
