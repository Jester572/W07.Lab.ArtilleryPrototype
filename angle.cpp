#include "angle.h"

/************************************
 * ANGLE
 ************************************/


 // default constructor
Angle::Angle() {
	this->radians = 0.0;
}

//Non-default constructor
Angle::Angle(double degrees) {
	this->radians = convertToRadians(degrees);
}

//copy constructor
Angle::Angle(const Angle& rhs) {
	this->radians = rhs.radians;
}

double Angle::getDegrees() const {
	double degrees = convertToDegrees(radians);
	return degrees;
};

double Angle::getRadians() const {
	double new_radians = normalize(radians);
	return new_radians;
};

void Angle::setDegrees(double degrees) {
	radians = normalize(convertToRadians(degrees));
};

void Angle::setRadians(double new_radians) {
	radians = normalize(new_radians);
};

void Angle::display(ostream& out) const {
	double degrees = getDegrees();
	out.setf(ios::fixed | ios::showpoint);
	out.precision(1);

	out << degrees << "degrees";
};

double Angle::normalize(double radian) const {
	// If it is a negative makes it a positive
	while (radian < 0)
	{
		radian +=  PI * 2;
	}
	// If Greater than 2 PI (360 degrees) sets it within 0, 2 PI
	while (radian >= PI * 2)
	{
		radian -= PI * 2;
	}

	return radian;
};

double Angle::convertToDegrees(double radians) const {
	double degrees = normalize(radians) * 360 / PI * 2;
	return degrees;
};

/*****************************************************************************
 * RADIANS_FROM_DEGREES
 * Convert from degrees to radians
 * INPUT:
 *   deg - Angle in degrees
 *****************************************************************************/
double Angle::radiansFromDegrees(const double& deg) {
	return (deg / 180.0) * PI;
};