//
// Created by jf1519 on 22/03/23.
//

#ifndef PIC_SEMI_IMPLICIT_VECTOR2_H
#define PIC_SEMI_IMPLICIT_VECTOR2_H

#include "Vector3.h"

template<typename T>
struct Vector2 {
    const unsigned int n;
    T* x; //x and y are contiguous in memory, so an increment of n leads to the next dimension.
    T* y;
    Vector2() = delete;
    explicit Vector2(unsigned int n);
    Vector2(const Vector2& other);
    Vector2(Vector2&& other) noexcept;
    ~Vector2();
    template <typename U> friend void swap(Vector2<U>& first, Vector2<U>& second);
    Vector2& operator=(Vector2 other);

    Vector2<T> operator+(const Vector2<T>& other) const;
    Vector2<T> operator+(const Vector3<T>& other) const; //Third column is ignored
    Vector2<T> operator-(const Vector2<T>& other) const;
    Vector2<T> operator-(const Vector3<T>& other) const; //Third column is ignored
    Vector2<T> operator*(double c) const;
    Vector2<T>& operator+=(const Vector2<T>& other);
    Vector2<T>& operator+=(const Vector3<T>& other); //Third column is ignored

    T* operator[](const unsigned int& dim); //Returns pointer to the numbered dimension
    const T* operator[](const unsigned int& dim) const; //Read-only version
};


#endif //PIC_SEMI_IMPLICIT_VECTOR2_H
