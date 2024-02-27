#include "drag.h"
#include <cmath>
#define PI atan(1) * 4  // PI

Drag::Drag()
{
	this->dragCoefficient = 0.3;
	this->dragDensity = 0.6;
}
Drag::Drag(double coefficient, double density)
{
	this->dragCoefficient = coefficient;
	this->dragDensity = density;
}

double Drag::getCoefficient()
{
	return dragCoefficient;
}

void Drag::setCoefficient(double newCoefficient)
{
	this->dragCoefficient = newCoefficient;
}

double Drag::getDensity()
{
	return dragDensity;
}
void Drag::setDensity(double newDensity)
{
	this->dragDensity = newDensity;
}

double Drag::computeDrag(double)
{
	const double radius = diameter / 2;

	double area = PI * pow(radius, 2);

	// force in newtons
	double drag = .5 * dragCoefficient * density * pow(velocity, 2) * area;

	return drag;
}