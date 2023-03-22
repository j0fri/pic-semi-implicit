//
// Created by jf1519 on 22/03/23.
//

#ifndef PIC_SEMI_IMPLICIT_VECTOR2_H
#define PIC_SEMI_IMPLICIT_VECTOR2_H


template<typename T>
struct Vector2 {
    const unsigned int n;
    T* x;
    T* y;
    Vector2() = delete;
    explicit Vector2(unsigned int n);
    Vector2(const Vector2& other);
    Vector2(Vector2&& other) noexcept;
    ~Vector2();
    template <typename U> friend void swap(Vector2<U>& first, Vector2<U>& second);
    Vector2& operator=(Vector2 other);

    Vector2<T> operator+(const Vector2<T>& other) const;
    Vector2<T> operator-(const Vector2<T>& other) const;
    Vector2<T> operator*(double c) const;
    Vector2<T>& operator+=(const Vector2<T>& other);
};


#endif //PIC_SEMI_IMPLICIT_VECTOR2_H
