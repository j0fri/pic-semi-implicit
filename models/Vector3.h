//
// Created by jf1519 on 22/03/23.
//

#ifndef PIC_SEMI_IMPLICIT_VECTOR3_H
#define PIC_SEMI_IMPLICIT_VECTOR3_H


template<typename T>
struct Vector3 {
    const unsigned int n;
    T* x;
    T* y;
    T* z;
    Vector3() = delete;
    explicit Vector3(unsigned int n);
    Vector3(const Vector3& other);
    Vector3(Vector3&& other) noexcept;
    ~Vector3();
    template <typename U> friend void swap(Vector3<U>& first, Vector3<U>& second);
    Vector3& operator=(Vector3 other);

    Vector3<T> operator+(const Vector3<T>& other) const;
    Vector3<T> operator-(const Vector3<T>& other) const;
    Vector3<T> operator*(double c) const;
    Vector3<T>& operator+=(const Vector3<T>& other);

    T* operator[](const unsigned int& dim); //Returns pointer to the numbered dimension
    const T* operator[](const unsigned int& dim) const; //Read-only version
};


#endif //PIC_SEMI_IMPLICIT_VECTOR3_H
