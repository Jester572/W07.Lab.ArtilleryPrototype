#include "velocity.h"

Velocity::Velocity() : dx(0.0), dy(0.0), mach(0.0)
{
}

Velocity::Velocity(double dx, double dy)
	: mach(0)
{
	this->dx = dx;
	this->dy = dy;
}

//Getters and setters
double Velocity::getDX()
{
	return dx;
}

void Velocity::setDX(double newDX)
{
	this->dx = newDX;
}

double Velocity::getDY()
{
	return dy;
}

void Velocity::setDY(double newDY)
{
	this->dx = newDY;
}

double Velocity::computeVelocity(double velocity, double acceleration, double timeInterval)
{
	// Find the new velocity.
	double newVelocity;
	newVelocity = velocity + (acceleration * timeInterval);
	return newVelocity;
}

double Velocity::updateVelocity(Acceleration& acceleration, double timeInterval, double altitude)
{

}