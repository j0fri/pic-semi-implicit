#include "Field2D3VConst.h"

template <typename T>
Field2D3VConst<T>::Field2D3VConst(const typename Config<T,2,3>::FieldConfig& fieldConfig,
                                  const typename Config<T,2,3>::BCConfig& bcConfig):
                                  Field2D3V<T>(fieldConfig, bcConfig){}

template <typename T>
void Field2D3VConst<T>::solveAndAdvance(T dt){
    std::copy(this->field, this->field+6*this->Ng, this->fieldT);
}

template class Field2D3VConst<float>;
template class Field2D3VConst<double>;