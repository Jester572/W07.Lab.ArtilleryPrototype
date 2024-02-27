#include <cmath>

/*****************************************************************************
 * GET_HORIZONTAL_SPEED
 * Compute horizontal component of speed from the angle in radians and the
 * overall speed.
 * INPUT:
 *   speed - Overall speed
 *   rad - Angle in radians
 *****************************************************************************/
double getHorizontalSpeed(const double& speed, const double& rad)
{
	return speed * sin(rad);
}

/*****************************************************************************
 * GET_VERTICAL_SPEED
 * INPUT:
 *   speed - Overall speed
 *   rad - Angle in radians
 *****************************************************************************/
double getVerticalSpeed(const double& speed, const double& rad)
{
	return speed * cos(rad);
}

/*****************************************************************************
 * GET_TOTAL_SPEED
 * INPUT:
 *   dx - Horizontal speed
 *   dy - Vertical speed
 *****************************************************************************/
double getTotalSpeed(const double& dx, const double& dy)
{
	return sqrt(dx * dx + dy * dy);
}