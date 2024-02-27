#include "gravity.h"

Gravity::Gravity()
{
	this->gravity = 
    {
   {0, -9.807},
   {1000, -9.804},
   {2000, -9.801},
   {3000, -9.797},
   {4000, -9.794},
   {5000, -9.791},
   {6000, -9.788},
   {7000, -9.785},
   {8000, -9.782},
   {9000, -9.779},
   {10000, -9.776},
   {15000, -9.761},
   {20000, -9.745},
   {25000, -9.730}
    };
}

/*****************************************************************************
 * LINEAR INTERPOLATION
 * curve fits between two points
 * INPUT:
 *   x0 - x value before the point interested
 *   y0 - y value before the point interested
 *   x1 - x value after the point interested
 *   y1 - y value after the point interrested
 *   x - current x value
 *****************************************************************************/
double linearInterpolation(double x0, double y0, double x1, double y1, double x)
{
    double y = y0 + ((x - x0) * ((y1 - y0) / (x1 - x0)));

    return y;
}

/*****************************************************************************
 * COMPUTE GRAVITY
 * finds gravity force
 * INPUT:
 *   altittude -  the current altitude
 *****************************************************************************/
double Gravity::computeGravity(int altitude)
{
    for (int i = 0; i <= 9000; i += 1000)
    {
        if (altitude >= i)
        {
            if (altitude < i + 1000)
            {
                return linearInterpolation(i, this->gravity.find(i)->second, i + 1000, this->gravity.find(i + 1000)->second, altitude);
            }
        }
    }
    for (int i = 10000; i < 25000; i += 5000)
    {
        if (altitude > i)
        {
            if (altitude <= i + 5000)
            {
                return linearInterpolation(i, this->gravity.find(i)->second, i + 5000, this->gravity.find(i + 5000)->second, altitude);
            }
        }
    }
    return -9.804;
}