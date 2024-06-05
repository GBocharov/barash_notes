#include "Quaternion.h"
#include <iostream>

Quaternion Quaternion::operator*(Quaternion A) const {
    Quaternion result{};
    result.r = r * A.r - i * A.i - j * A.j - k * A.k;
    result.i = r * A.i + i * A.r + j * A.k - k * A.j;
    result.j = r * A.j - i * A.k + j * A.r + k * A.i;
    result.k = r * A.k + i * A.j - j * A.i + k * A.r;
    return result;
}

Matrix Quaternion::toMatrix() const {
    Matrix result{};
    result.values[0][0] = 1 - 2 * j * j - 2 * k * k;
    result.values[0][1] = 2 * i * j - 2 * r * k;
    result.values[0][2] = 2 * i * k + 2 * r * j;
    result.values[1][0] = 2 * i * j + 2 * r * k;
    result.values[1][1] = 1 - 2 * i * i - 2 * k * k;
    result.values[1][2] = 2 * j * k - 2 * r * i;
    result.values[2][0] = 2 * i * k - 2 * r * j;
    result.values[2][1] = 2 * j * k + 2 * r * i;
    result.values[2][2] = 1 - 2 * i * i - 2 * j * j;
    return result;
}

Quaternion Quaternion::normalize() const {
    double length = sqrt(r * r + i * i + j * j + k * k);
    if (length != 0) {
        return Quaternion{ r / length, i / length, j / length, k / length };
    }
    return Quaternion{};
}

Quaternion Quaternion::operator*(double h) const {
    Quaternion result{};
    result.r = r * h;
    result.i = i * h;
    result.j = j * h;
    result.k = k * h;
    return result;
}

Quaternion Quaternion::operator+(Quaternion other) const {

    Quaternion result{};
    result.r = r + other.r;
    result.i = i + other.i;
    result.j = j + other.j;
    result.k = k + other.k;
    return result;

}