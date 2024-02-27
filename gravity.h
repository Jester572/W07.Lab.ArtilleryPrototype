#pragma once
#include <map>
class Gravity
{
private:
	std::map<int, double> gravity;

public:
	double computeGravity(int altitude);

};

