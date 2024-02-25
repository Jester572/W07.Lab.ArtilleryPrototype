/*************************************************************
 * 1. Name:
 *      Nathan Bird/Zach Wilson/Jesse Earley
 * 2. Assignment Name:
 *      Lab 07: Artillery Prototype
 * 3. Assignment Description:
 *      Compute the distance and hang time of an M795 posProjectile
 *      when fired from the M777 howitzer.
 * 4. What was the hardest part? Be as specific as possible.
 *      ??
 * 5. How long did it take for you to complete the assignment?
 *      1.5h + ?? + ??
 *****************************************************************/

#include <cassert>      // for ASSERT
#include <cmath>       // for math functions
#include "position.h"   // for POSITION
#include <map>
using namespace std;

#define PI atan(1) * 4  // PI

/*****************************************************************************
 * RADIANS_FROM_DEGREES
 * Convert from degrees to radians
 * INPUT:
 *   deg - Angle in degrees
 *****************************************************************************/
double radiansFromDegrees(const double & deg)
{
   return (deg / 180.0) * PI;
}

/*****************************************************************************
 * GET_HORIZONTAL_SPEED
 * Compute horizontal component of speed from the angle in radians and the
 * overall speed.
 * INPUT:
 *   speed - Overall speed
 *   rad - Angle in radians
 *****************************************************************************/
double getHorizontalSpeed(const double & speed, const double & rad)
{
   return speed * sin(rad);
}

/*****************************************************************************
 * GET_VERTICAL_SPEED
 * INPUT:
 *   speed - Overall speed
 *   rad - Angle in radians
 *****************************************************************************/
double getVerticalSpeed(const double & speed, const double & rad)
{
   return speed * cos(rad);
}

/*****************************************************************************
 * GET_TOTAL_SPEED
 * INPUT:
 *   dx - Horizontal speed
 *   dy - Vertical speed
 *****************************************************************************/
double getTotalSpeed(const double & dx, const double & dy)
{
   return sqrt(dx * dx + dy * dy);
}

/*****************************************************************************
 * COMPUTE_VELOCITY
 * Gets the new velocity from old velocity, acceleration, and time
 * INPUT:
 *   vel - Old velocity
 *   accel - Acceleration
 *   time - Time
 *****************************************************************************/
double computeVelocity(double vel, double accel, double time)
{
   return vel + (accel * time);
}

/*****************************************************************************
 * COMPUTE_DISTANCE
 * Gets the new distance from old distance, velocity, acceleration, and time
 * INPUT:
 *   s1 - Old distance
 *   vel - Velocity
 *   accel - Acceleration
 *   time - Time
 *****************************************************************************/
double computeDistance(double s1, double vel, double accel, double time)
{
   return s1 + (vel * time) + (0.5 * (accel * (time * time)));
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
 *   gravity_map - A mapped key, value pairs of altitude and gravity forces
 *   altittude -  the current altitude
 *****************************************************************************/
double computeGravity(map<int, double> gravity_map, double altitude)
{
    for (int i = 0; i <= 9000; i += 1000)
    {
        if (altitude >= i) 
        {
            if (altitude < i + 1000)
            {
                return linearInterpolation(i, gravity_map.find(i)->second, i + 1000, gravity_map.find(i + 1000)->second, altitude);
            }
        }
    }
    for (int i = 10000; i < 25000; i += 5000)
    {
        if (altitude > i)
        {
            if (altitude <= i + 5000)
            {
                return linearInterpolation(i, gravity_map.find(i)->second, i + 5000, gravity_map.find(i + 5000)->second, altitude);
            }
        }
    }
    return -9.804;
}


double computeDrag(double velocity, double diameter) {
    const double dragCoefficient = 0.3;
    const double density = 0.6;
    //diameter in meters

    const double radius = diameter / 2;

    double area = PI * pow(radius, 2);

    // force in newtons
    double drag = .5 * dragCoefficient * density * pow(velocity, 2) * area;

    return drag;
}

/*****************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
   Position posProjectile(0.0, 0.0);  // starting distance, altitude (x, y)

   double angleDeg = 75.0;  // degrees
   double angleRad = radiansFromDegrees(angleDeg);  // radians

   const double muzzleVel = 827.0;  // m/s
   const double projectileWeight = 46.7;  // kg
   const double projectileDiameter = 0.15489; // m

   // Acceleration
   double ddx;
   double ddy;
   std::map<int, double> gravity = {
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

   // Break down muzzle velocity into horizontal and vertical, applying angle
   double dx = getHorizontalSpeed(muzzleVel, angleRad);
   double dy = getVerticalSpeed(muzzleVel, angleRad);

   double dt = 0.01;
   double x;
   double y;

   // Total initial velocity
   double initialVelocity = getTotalSpeed(dx, dy);



   ///* 1. Inertia: No drag, no accel, no gravity; so, dx and dy won't change */
   //for (int t = 0; t < 20; t++)
   //{
   //   // Update position based on constant speed
   //   posProjectile.addMetersX(dx);
   //   posProjectile.addMetersY(dy);
   //}
   //cout << "Distance: " << posProjectile.getMetersX() << "m  Altitude: " << posProjectile.getMetersY() << "m\n";



   /* 2. Acceleration: introduce constant gravity */
   double totalTime = 0;
   while(posProjectile.getMetersY() >= 0)
   {
      // all the different accelerations
      double computedGravity = computeGravity(gravity, posProjectile.getMetersY());
      double dragx = -1 * computeDrag(dx, projectileDiameter);
      double dragy = -1 * computeDrag(dy, projectileDiameter);

      // acceleration split between x and y
      ddx = dragx;
      ddy = computedGravity + dragy;

      // x and y components of velocity
      dx = computeVelocity(dx, ddx, dt);
      dy = computeVelocity(dy, ddy, dt);

      // positions
      x = computeDistance(posProjectile.getMetersX(), dx, ddx, dt);
      y = computeDistance(posProjectile.getMetersY(), dy, ddy, dt);
      // Update position
      posProjectile.setMetersX(x);
      posProjectile.setMetersY(y);
      totalTime += dt;
   }
   cout << "Distance: " << posProjectile.getMetersX() << "m  Altitude: " << posProjectile.getMetersY() << "m Hang Time: " << totalTime << "s" << endl;

   return EXIT_SUCCESS;
}
