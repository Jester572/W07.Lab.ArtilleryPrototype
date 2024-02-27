#pragma once
#include <math.h>    // for floor()
#include <iostream>  // for cout
#include <cassert>   // for assert()
using namespace std;

#define PI atan(1) * 4  // PI

class Angle
{
private:
	double radians; //angle

	double normalize(double angle) const;

	double convertToDegrees(double radians) const;

	double convertToRadians(double degrees) const;

public:
	// default constructor
	Angle();

	//Non-default constructor
	Angle(double radians);

	//copy constructor
	Angle(const Angle& rhs);

	double getDegrees() const;

	double getRadians() const;

	void setDegrees(double degrees);

	void setRadians(double new_radians);

	void display(ostream& out) const;

};

