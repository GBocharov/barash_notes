#ifndef VECTOR_H
#define VECTOR_H
#include "Matrix.h"

struct Matrix;

struct Vector {
    double x;
    double y;
    double z;


    Vector operator*(const Vector& other) const;
    Vector operator*(const Matrix& matrix) const;
    double dot(const Vector& other) const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(const double d) const;
    Vector operator/(const double d) const;

    double abs();

    void print();
};


#endif //VECTOR_H