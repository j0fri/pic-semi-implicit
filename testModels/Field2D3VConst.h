#ifndef PIC_SEMI_IMPLICIT_FIELD2D3VCONST_H
#define PIC_SEMI_IMPLICIT_FIELD2D3VCONST_H

#include "../models/Field2D3V.h"
#include "../models/Config.h"

//2D3V field which is constant after initialisation for debug purposes.
//It runs all accumulation operations but not system solve.
template <typename T>
class Field2D3VConst: public Field2D3V<T> {
public:
    Field2D3VConst() = delete;
    explicit Field2D3VConst(const typename Config<T,2,3>::FieldConfig& fieldConfig,
                            const typename Config<T,2,3>::BCConfig& bcConfig);
    Field2D3VConst(const Field2D3VConst& other) = delete;
    ~Field2D3VConst() = default;
private:
    void solveAndAdvance(T dt);
};

#endif //PIC_SEMI_IMPLICIT_FIELD2D3VCONST_H