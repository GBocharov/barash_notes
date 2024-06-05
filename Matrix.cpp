#include "Matrix.h"
#include <iostream>

Matrix Matrix::operator*(const Matrix& other) const {
    Matrix result{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result.values[i][j] = 0;
            for (int k = 0; k < 3; ++k) result.values[i][j] += values[i][k] * other.values[k][j];
        }
    }
    return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
    Matrix result{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result.values[i][j] = values[i][j] + other.values[i][j];
        }
    }
    return result;
}


Vector Matrix::operator*(const Vector& vector) {


    return Vector{
            values[0][0] * vector.x + values[0][1] * vector.y + values[0][2] * vector.z,
            values[1][0] * vector.x + values[1][1] * vector.y + values[1][2] * vector.z,
            values[2][0] * vector.x + values[2][1] * vector.y + values[2][2] * vector.z
    };
}

Matrix Matrix::operator/(const double d) const {
    Matrix result{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result.values[i][j] = values[i][j] /d;
        }
    }
    return result;
}



Matrix Matrix::transpose() {
    Matrix result{0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            result.values[i][j] = values[j][i];
    return result;
}

void Matrix::printMatr()
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
        {
            std::cout << values[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


Matrix Matrix::Inverse()
{
    Matrix res;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int i1, i2, j1, j2;
            if (i == 0) {
                i1 = 1;  i2 = 2;
            }
            else if (i == 1) {
                i1 = 0;  i2 = 2;
            }
            else {
                i1 = 0;  i2 = 1;
            }
            if (j == 0) {
                j1 = 1;  j2 = 2;
            }
            else if (j == 1) {
                j1 = 0;  j2 = 2;
            }
            else {
                j1 = 0;  j2 = 1;
            }
            if ((i + j) % 2) {
                if ((values[i1][j1] * values[i2][j2] - values[i1][j2] * values[i2][j1] != 0))
                    res.values[j][i] = -(values[i1][j1] * values[i2][j2] - values[i1][j2] * values[i2][j1]);
                else
                    res.values[j][i] = 0;
            }
            else {
                res.values[j][i] = (values[i1][j1] * values[i2][j2] - values[i1][j2] * values[i2][j1]);
            }
        }
    }

    return (res/det());
}


double Matrix::det() {
    double det = 0.0;

    for (int i = 0; i < 3; i++) {
        int i1, i2, j1, j2;
        if (i == 0) {
            i1 = 1;  i2 = 2;
        }
        else if (i == 1) {
            i1 = 0;  i2 = 2;
        }
        else {
            i1 = 0;  i2 = 1;
        }
        j1 = 1; j2 = 2;
        if (i % 2) {
            det += -(values[i1][j1] * values[i2][j2] - values[i1][j2] * values[i2][j1]) * values[i][0];
        }
        else {
            det += (values[i1][j1] * values[i2][j2] - values[i1][j2] * values[i2][j1]) * values[i][0];
        }
    }

    return det;
}
