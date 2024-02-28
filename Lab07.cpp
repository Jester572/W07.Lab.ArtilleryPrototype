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

/******************************************
 * NORMALIZE
 *    returns the "unwrapped" value of rad
 *****************************************/
double normalize(double rad)
{
   double norm_rad = rad;

   if (norm_rad > 0)
   {
      while (norm_rad > PI * 2.0)
         norm_rad -= PI * 2;
   }
   else
   {
      while (norm_rad < 0)
         norm_rad += PI * 2;
   }

   return norm_rad;
}

/*****************************************************************************
 * GET_HORIZONTAL_COMPONENT
 * Compute horizontal component of total from the angle in radians and the
 * overall total.
 * INPUT:
 *   total - Overall total
 *   rad - Angle in radians
 *****************************************************************************/
double getHorizontalComponent(const double & total, const double & rad)
{
   return total * sin(rad);
}

/*****************************************************************************
 * GET_VERTICAL_COMPONENT
 * INPUT:
 *   total - Overall total
 *   rad - Angle in radians
 *****************************************************************************/
double getVerticalComponent(const double & total, const double & rad)
{
   return total * cos(rad);
}

/*****************************************************************************
 * GET_TOTAL_COMPONENT
 * INPUT:
 *   x - Horizontal total
 *   y - Vertical total
 *****************************************************************************/
double getTotalComponent(const double & x, const double & y)
{
   return sqrt(x * x + y * y);
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
 *   key0 - x value before the point interested
 *   value0 - y value before the point interested
 *   key1 - x value after the point interested
 *   value1 - y value after the point interrested
 *   x - current x value
 *****************************************************************************/
double linearInterpolation(double x0, double y0, double x1, double y1, double x)
{
    double y = y0 + ((x - x0) * ((y1 - y0) / (x1 - x0)));
    
    return y;
}

/*****************************************************************************
 * INTERPOLATE_VALUE_FROM_MAP
 * Finds in the inputs (keys) of the map where inputValue would be and uses the
 * two closest points (key, value) to interpolate its outputValue.
 * INPUT:
 *   tableMap - A mapped key, value pairs
 *   inputValue - the input for the map
 *****************************************************************************/
double interpolateValueFromMap(map<double, double> tableMap, double inputValue)
{
   double key0;
   double value0;
   double key1;
   double value1;
   for (auto it = tableMap.begin(); it != tableMap.end(); it)
   {
      key0 = it->first;
      value0 = it->second;

      if (++it == tableMap.end())
         return it->second;  // Upper bound

      key1 = it->first;
      value1 = it->second;

      if (inputValue >= key0 && inputValue < key1)
         return linearInterpolation(key0, value0, key1, value1, inputValue);
   }

   return tableMap.begin()->second;  // Lower bound
}

double computeDrag(const double& velocity, const double& dragCoefficient, const double& density, const double& diameter) {
    //const double dragCoefficient = 0.3;
    //diameter in meters

    const double radius = diameter / 2;

    double area = PI * pow(radius, 2);

    // force in newtons
    double drag = .5 * dragCoefficient * density * pow(velocity, 2) * area;

    return drag;
}

double computeAcceleration(const double& force, const double& mass)
{
   return force / mass;
}

/*****************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
   Position posProjectile(0.0, 0.0);  // starting distance, inputValue (x, y)

   double angleDeg = 75.0;  // degrees
   double angleRad = radiansFromDegrees(angleDeg);  // radians

   const double muzzleVel = 827.0;  // m/s
   const double projectileWeight = 46.7;  // kg
   const double projectileDiameter = 0.15489; // m

   // Acceleration
   double ddx;
   double ddy;

   // Altitude to gravity (m/s^2)
   std::map<double, double> gravityMap = {
       {0.0, -9.807},
       {1000.0, -9.804},
       {2000.0, -9.801},
       {3000.0, -9.797},
       {4000.0, -9.794},
       {5000.0, -9.791},
       {6000.0, -9.788},
       {7000.0, -9.785},
       {8000.0, -9.782},
       {9000.0, -9.779},
       {10000.0, -9.776},
       {15000.0, -9.761},
       {20000.0, -9.745},
       {25000.0, -9.730}
   };

   // Altitude to air density (kg/m^2)
   map<double, double> airDensityMap = {
      {0.0, 1.2250000},
      {1000.0, 1.1120000},
      {2000.0, 1.0070000},
      {3000.0, 0.9093000},
      {4000.0, 0.8194000},
      {5000.0, 0.7364000},
      {6000.0, 0.6601000},
      {7000.0, 0.5900000},
      {8000.0, 0.5258000},
      {9000.0, 0.4671000},
      {10000.0, 0.4135000},
      {15000.0, 0.1948000},
      {20000.0, 0.0889100},
      {25000.0, 0.0400800},
      {30000.0, 0.0184100},
      {40000.0, 0.0039960},
      {50000.0, 0.0010270},
      {60000.0, 0.0003097},
      {70000.0, 0.0000828},
      {80000.0, 0.0000185}
   };

   // Altitude to speed of sound (mach) (m/s)
   map<double, double> machMap = {
      {0.0, 340.0},
      {1000.0, 336.0},
      {2000.0, 332.0},
      {3000.0, 328.0},
      {4000.0, 324.0},
      {5000.0, 320.0},
      {6000.0, 316.0},
      {7000.0, 312.0},
      {8000.0, 308.0},
      {9000.0, 303.0},
      {10000.0, 299.0},
      {15000.0, 295.0},
      {20000.0, 295.0},
      {25000.0, 295.0},
      {30000.0, 305.0},
      {40000.0, 324.0}
   };

   // Mach to drag coefficient (force component)
   map<double, double> dragCoefficientMap = {
      {0.300, 0.1629},
      {0.500, 0.1659},
      {0.700, 0.2031},
      {0.890, 0.2597},
      {0.920, 0.3010},
      {0.960, 0.3287},
      {0.980, 0.4002},
      {1.000, 0.4258},
      {1.020, 0.4335},
      {1.060, 0.4483},
      {1.240, 0.4064},
      {1.530, 0.3663},
      {1.990, 0.2897},
      {2.870, 0.2297},
      {2.890, 0.2306},
      {5.000, 0.2656}
   };

   // Break down muzzle velocity into horizontal and vertical, applying angle
   double dx = getHorizontalComponent(muzzleVel, angleRad);
   double dy = getVerticalComponent(muzzleVel, angleRad);

   double dt = 0.01;
   double x;
   double y;

   // Total initial velocity
   double initialVelocity = getTotalComponent(dx, dy);

   double totalVel;
   double mach;
   double mVel;

   double airDensity;
   double dragCoefficient;
   double dragForce;
   double acc;
   double computedGravity;

   double theta;
   double ddyDrag;

   double totalTime = 0;
   while(posProjectile.getMetersY() >= 0)
   {
       // totals
       totalVel = getTotalComponent(dx, dy);
       mach = interpolateValueFromMap(machMap, posProjectile.getMetersY());  // Speed of sound at current altitude
       mVel = totalVel / mach;  // Speed of projectile in relation to speed of sound

       airDensity = interpolateValueFromMap(airDensityMap, posProjectile.getMetersY());
       dragCoefficient = interpolateValueFromMap(dragCoefficientMap, mVel);
       dragForce = computeDrag(totalVel, dragCoefficient, airDensity, projectileDiameter);
       acc = computeAcceleration(dragForce, projectileWeight);
       computedGravity = interpolateValueFromMap(gravityMap, posProjectile.getMetersY());

       // drag (acceleration) components
       theta = normalize(atan2(dx, dy) + PI);  // adding PI radians gets the opposite of the computed angle
       ddx = getHorizontalComponent(acc, theta);
       ddyDrag = getVerticalComponent(acc, theta);
       ddy = computedGravity + ddyDrag;

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
