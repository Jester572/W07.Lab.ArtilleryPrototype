#pragma once
#include "forces.h"
class Acceleration
{
private:
	double ddx;
	double ddy;

public:
	Acceleration();

	// Methods
	void   updateAcceleration(Force& force, double mass);
	double computeAcceleration(double force, double mass);
	
	// Getters and Setters
	double getDDX();
	double getDDY();

	void setDDX();
	void setDDY();
};

