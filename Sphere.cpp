#include "Sphere.h"
#include <iostream>


Sphere::Sphere(Vector position, double m, double r) {
	mass = m;
	coord = position;
	Radius = r;
	I = {0,0,0,0,0,0,0,0,0};
	//std::cout << "R = " << R << "mass = " << mass << std::endl;
	I.values[0][0] = (2.0 / 5.0) * mass * Radius * Radius; 
	I.values[1][1] = (2.0 / 5.0) * mass * Radius * Radius;
	I.values[2][2] = (2.0 / 5.0) * mass * Radius * Radius;
}

void Sphere::Compute(Vector v) // Teorema  G - Shteiner Translate tensor
{
	double a = coord.x - v.x; // 
	double b = coord.y - v.y;
	double c = coord.z - v.z;

	I.values[0][0] += mass * (b * b + c * c);
	I.values[1][1] += mass * (a * a + c * c);
	I.values[2][2] += mass * (a * a + b * b);

	I.values[0][1] -= mass * a * b;
	I.values[0][2] -= mass * a * c;
	I.values[1][2] -= mass * b * c;

	I.values[1][0] -= mass * a * b;
	I.values[2][0] -= mass * a * c;
	I.values[2][1] -= mass * b * c;

}