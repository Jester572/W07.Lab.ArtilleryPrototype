#pragma once

#include "acceleration.h"
#include <iostream>
using namespace std;

class Velocity
{
private:
	double dx;
	double dy;
	double mach;


public:
	//Default Constructor
	Velocity();
	
	//Parameterized Constructor
	Velocity(double dx, double dy);

	void addX(double accel, double time);
	void addY(double accel, double time);

	double computeVelocity(double vel, double accel, double time);
	void updateVelocity(Acceleration &acceleration, double timeInterval, double altitude);
	double	computeMach(double altitude);

	//Getters and Setters
	double  getDX();
	double getDY();
	double getVelocity();
	double getMach();

	void setDX(double newDX);
	void setDY(double newDY);

};

