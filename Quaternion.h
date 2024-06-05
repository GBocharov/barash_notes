#ifndef WITCH2_QUATERNION_H
#define WITCH2_QUATERNION_H

#include "Matrix.h"
#include "cmath"

struct Quaternion {
    double r;
    double i;
    double j;
    double k;

    Quaternion operator*(Quaternion A) const;
    Quaternion operator*(double h) const;
    Quaternion operator+(Quaternion other) const;
    Matrix toMatrix() const;
    Quaternion normalize() const;
};


#endif //WITCH2_QUATERNION_H