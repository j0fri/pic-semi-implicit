//
// Created by jf1519 on 22/03/23.
//

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "Vector2.h"


template<typename T>
Vector2<T>::Vector2(unsigned int n): n{n} {
    try{
        x = new T[2*n];
        y = x + n;
    }catch(const std::bad_alloc& e){
        std::cerr << "Exception thrown in memory allocation in Vector2" << std::endl;
        throw;
    }
}

template<typename T>
Vector2<T>::Vector2(const Vector2 &other): n{other.n} {
    x = new T[n];
    y = new T[n];
    std::copy(other.x, other.x+n, x);
    std::copy(other.y, other.y+n, y);
}

template<typename T>
Vector2<T>::Vector2(Vector2 &&other) noexcept : n{other.n}{
    this->x = other.x;
    this->y = other.y;
}

template<typename T>
Vector2<T>::~Vector2() {
    delete[] x;
}

template<typename T>
Vector2<T> &Vector2<T>::operator=(Vector2 other) {
    swap(*this,other);
    return *this;
}

template<typename U>
void swap(Vector2<U>& first, Vector2<U>& second) {
    using std::swap;
    if(first.n != second.n){
        throw std::runtime_error("Can only swap vectors with the same length.");
    }
    std::swap(first.x, second.x);
    std::swap(first.y, second.y);
}
template void swap<float>(Vector2<float>& first, Vector2<float>& second);
template void swap<double>(Vector2<double>& first, Vector2<double>& second);
template void swap<int>(Vector2<int>& first, Vector2<int>& second);
template void swap<unsigned int>(Vector2<unsigned int>& first, Vector2<unsigned int>& second);


template<typename T>
Vector2<T> Vector2<T>::operator+(const Vector2 &other) const {
    if(this->n != other.n){
        throw std::runtime_error("Can only operate on vectors with the same length.");
    }
    Vector2<T> out(n);
    for(unsigned int i = 0; i < n; ++i){
        out.x[i] = this->x[i] + other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        out.y[i] = this->y[i] + other.y[i];
    }
    return out;
}

template<typename T>
Vector2<T> Vector2<T>::operator+(const Vector3<T> &other) const {
    if(this->n != other.n){
        throw std::runtime_error("Can only operate on vectors with the same length.");
    }
    Vector2<T> out(n);
    for(unsigned int i = 0; i < n; ++i){
        out.x[i] = this->x[i] + other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        out.y[i] = this->y[i] + other.y[i];
    }
    return out;
}

template<typename T>
Vector2<T> Vector2<T>::operator-(const Vector2 &other) const {
    if(this->n != other.n){
        throw std::runtime_error("Can only operate on vectors with the same length.");
    }
    Vector2<T> out(n);
    for(unsigned int i = 0; i < n; ++i){
        out.x[i] = this->x[i] - other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        out.y[i] = this->y[i] - other.y[i];
    }
    return out;
}

template<typename T>
Vector2<T> Vector2<T>::operator-(const Vector3<T> &other) const {
    if(this->n != other.n){
        throw std::runtime_error("Can only operate on vectors with the same length.");
    }
    Vector2<T> out(n);
    for(unsigned int i = 0; i < n; ++i){
        out.x[i] = this->x[i] - other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        out.y[i] = this->y[i] - other.y[i];
    }
    return out;
}

template<typename T>
Vector2<T> Vector2<T>::operator*(double c) const{
    Vector2<T> out(n);
    for(unsigned int i = 0; i < n; ++i){
        out.x[i] = c*x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        out.y[i] = c*y[i];
    }
    return out;
}

template<typename T>
Vector2<T> &Vector2<T>::operator+=(const Vector2<T> &other) {
    for(unsigned int i = 0; i < n; ++i){
        x[i] += other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        y[i] += other.y[i];
    }
    return *this;
}

template<typename T>
Vector2<T> &Vector2<T>::operator+=(const Vector3<T> &other) {
    for(unsigned int i = 0; i < n; ++i){
        x[i] += other.x[i];
    }
    for(unsigned int i = 0; i < n; ++i){
        y[i] += other.y[i];
    }
    return *this;
}

template<typename T>
T* Vector2<T>::operator[](const unsigned int &dim){
    if(dim == 0){
        return x;
    }else if(dim == 1){
        return y;
    }
    throw std::invalid_argument("Vector 2 only has dimensions 0 and 1");
}

template<typename T>
const T *Vector2<T>::operator[](const unsigned int &dim) const {
    if(dim == 0){
        return static_cast<const T*>(x);
    }else if(dim == 1){
        return static_cast<const T*>(y);
    }
    throw std::invalid_argument("Vector 2 only has dimensions 0 and 1");
}

template struct Vector2<float>;
template struct Vector2<double>;
template struct Vector2<int>;
template struct Vector2<unsigned int>;