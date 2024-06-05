
#pragma once

#include "Matrix.h"
#include "Vector.h"

struct Sphere
{
	double mass;
	double Radius; //rad
	Vector coord; // coord
	Matrix I; // tensor
	Sphere(Vector position, double m, double r);

	void Compute(Vector v); // // Teorema  G - Shteiner Translate tensor

};

