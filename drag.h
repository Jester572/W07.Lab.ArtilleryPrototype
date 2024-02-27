#pragma once
class Drag
{
private:
	double dragCoefficient;
	double dragDensity;

public:
	Drag();
	Drag(double coefficient, double density);

	double getCoefficient();
	void setCoefficient(double newCoefficient);

	double getDensity();
	void setDensity(double newDensity);

	double computeDrag(double velocity, double diameter)

};

