#ifndef WITCHMATH_MATRIX_H
#define WITCHMATH_MATRIX_H

#include "Vector.h"

struct Vector;

struct Matrix {
    double values[3][3];
    Matrix operator*(const Matrix& other) const;
    Matrix operator/(const double d) const;
    Matrix operator+(const Matrix& other) const;
    Vector operator*(const Vector& vector);
    Matrix transpose();
    Matrix Inverse();
    void printMatr();
    double det();
};


#endif //WITCHMATH_MATRIX_H