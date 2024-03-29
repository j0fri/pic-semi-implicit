#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "Vector3.h"


template<typename T>
Vector3<T>::Vector3(unsigned int n): n{n} {
    try{
        x = new T[3*n];
        y = x + n;
        z = y + n;
    }catch(const std::bad_alloc& e){
        std::cerr << "Exception thrown in memory allocation in Vector3" << std::endl;
        throw;
    }
}

template<typename T>
Vector3<T>::Vector3(const Vector3 &other): n{other.n} {
    try{
        x = new T[3*n];
        y = x + n;
        z = y + n;
    }catch(const std::bad_alloc& e){
        std::cerr << "Exception thrown in memory allocation in Vector3" << std::endl;
        throw;
    }
    std::copy(other.x, other.x+n, x);
    std::copy(other.y, other.y+n, y);
    std::copy(other.z, other.z+n, z);
}

template<typename T>
Vector3<T>::Vector3(Vector3 &&other) noexcept : n{other.n}{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    other.x = nullptr;
    other.y = nullptr;
    other.z = nullptr;
}

template<typename T>
Vector3<T>::~Vector3() {
    delete[] x;
}

template<typename T>
Vector3<T> &Vector3<T>::operator=(Vector3 other) {
    swap(*this,other);
    return *this;
}

template<typename U>
void swap(Vector3<U>& first, Vector3<U>& second) {
    using std::swap;
    if(first.n != second.n){
        throw std::runtime_error("Can only swap vectors with the same length.");
    }
    std::swap(first.x, second.x);
    std::swap(first.y, second.y);
    std::swap(first.z, second.z);
}
template void swap<float>(Vector3<float>& first, Vector3<float>& second);
template void swap<double>(Vector3<double>& first, Vector3<double>& second);
template void swap<int>(Vector3<int>& first, Vector3<int>& second);
template void swap<unsigned int>(Vector3<unsigned int>& first, Vector3<unsigned int>& second);


template<typename T>
Vector3<T> Vector3<T>::operator+(const Vector3 &other) const {
    if(this->n != other.n){
        throw std::runtime_error("Can only operate on vectors with the same length.");
    }
    Vector3<T> out(n);
//    for(unsigned int i = 0; i < n; ++i){
//        out.x[i] = this->x[i] + other.x[i];
//    }
//    for(unsigned int i = 0; i < n; ++i){
//        out.y[i] = this->y[i] + other.y[i];
//    }
//    for(unsigned int i = 0; i < n; ++i){
//        out.z[i] = this->z[i] + other.z[i];
//    }
    auto thisPtr = x-1;
    auto otherPtr = other.x-1;
    auto outPtr = out.x-1;
    for(unsigned int i = 0; i < 3*n; ++i){
        *(++outPtr) = *(++thisPtr) + *(++otherPtr);
    }
    return out;
}

template<typename T>
Vector3<T> Vector3<T>::operator-(const Vector3 &other) const {
    if(this->n != other.n){
        throw std::runtime_error("Can only operate on vectors with the same length.");
    }
    Vector3<T> out(n);
    //Unoptimised code
//    for(unsigned int i = 0; i < n; ++i){
//        out.x[i] = this->x[i] - other.x[i];
//    }
//    for(unsigned int i = 0; i < n; ++i){
//        out.y[i] = this->y[i] - other.y[i];
//    }
//    for(unsigned int i = 0; i < n; ++i){
//        out.z[i] = this->z[i] - other.z[i];
//    }
    auto thisPtr = x-1;
    auto otherPtr = other.x-1;
    auto outPtr = out.x-1;
    const auto last = x+3*n;
    while(++thisPtr != last){
        *(++outPtr) = *thisPtr - *(++otherPtr);
    }
    return out;
}

template<typename T>
Vector3<T> Vector3<T>::operator*(T c) const{
    Vector3<T> out(n);
    //Unoptimised code:
//    for(unsigned int i = 0; i < n; ++i){
//        out.x[i] = C*x[i];
//    }
//    for(unsigned int i = 0; i < n; ++i){
//        out.y[i] = C*y[i];
//    }
//    for(unsigned int i = 0; i < n; ++i){
//        out.z[i] = C*z[i];
//    }
    auto thisPtr = x-1;
    auto outPtr = out.x-1;
    for(unsigned int i = 0; i < 3*n; ++i){
        *(++outPtr) = c*(*(++thisPtr));
    }
    return out;
}

template<typename T>
Vector3<T> &Vector3<T>::operator+=(const Vector3<T> &other) {
    for(unsigned int i = 0; i < n; ++i){
        x[i] += other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        y[i] += other.y[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        z[i] += other.z[i];
    }
    return *this;
}

template<typename T>
T* Vector3<T>::operator[](const unsigned int &dim){
    if(dim == 0){
        return x;
    }else if(dim == 1){
        return y;
    }else if(dim == 2){
        return z;
    }
    throw std::invalid_argument("Vector 3 only has dimensions 0, 1 and 2");
}

template<typename T>
const T *Vector3<T>::operator[](const unsigned int &dim) const {
    if(dim == 0){
        return static_cast<const T*>(x);
    }else if(dim == 1){
        return static_cast<const T*>(y);
    }else if(dim == 2){
        return static_cast<const T*>(z);
    }
    throw std::invalid_argument("Vector 3 only has dimensions 0, 1 and 2");
}

template struct Vector3<float>;
template struct Vector3<double>;
template struct Vector3<int>;
template struct Vector3<unsigned int>;
