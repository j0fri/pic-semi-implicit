//
// Created by jf1519 on 22/03/23.
//
#include <iostream>
#include "../models/Vector2.h"
#include "../models/Vector3.h"


template<typename T>
void printVector2(const Vector2<T>& v, unsigned int n){
    for(unsigned int i = 0; i < n; ++i){
        std::cout << v.x[i] << " " << v.y[i] << std::endl;
    }
}

template<typename T>
void printVector3(const Vector3<T>& v, unsigned int n){
    for(unsigned int i = 0; i < n; ++i){
        std::cout << v.x[i] << " " << v.y[i] << " " << v.z[i] << std::endl;
    }
}


int main(){
    {
        int n = 5;
        Vector2<double> v1(n);
        for (int i = 0; i < n; ++i) {
            v1.x[i] = i;
            v1.y[i] = n - i;
        }
        std::cout << "x: " << std::endl;
        printVector2(v1, n);
        std::cout << "y=x: " << std::endl;
        auto v2 = v1;
        printVector2(v2, n);
        std::cout << "x+y: " << std::endl;
        printVector2(v1 + v2, n);
        std::cout << "z = -2y: " << std::endl;
        auto v3 = v2 * (-2);
        printVector2(v3, n);
        v3 += v2;
        std::cout << "z+=y, z: " << std::endl;
        printVector2(v3, n);
    }
    {
        int n = 5;
        Vector3<double> v1(n);
        for (int i = 0; i < n; ++i){
            v1.x[i] = i;
            v1.y[i] = n-i;
            v1.z[i] = 36;
        }
        std::cout << "x: " << std::endl;
        printVector3(v1,n);
        std::cout << "y=x: " << std::endl;
        auto v2 = v1;
        printVector3(v2,n);
        std::cout << "x+y: " << std::endl;
        printVector3(v1+v2,n);
        std::cout << "z = -2y: " << std::endl;
        auto v3 = v2*(-2);
        printVector3(v3,n);
        v3 += v2;
        std::cout << "z+=y, z: " << std::endl;
        printVector3(v3,n);
    }
}