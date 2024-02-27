#pragma once
#include "position.h"
#include "angle.h"
#include "acceleration.h"
#include "velocity.h"
#include "gravity.h"
#include "drag.h"

class Projectile
{
private:
	double mass;
	double surfaceArea;
	Position position;
	Angle angle;
	Acceleration acceleration;
	Velocity velocity;
	Gravity gravity;
	Drag drag;

public:


};

