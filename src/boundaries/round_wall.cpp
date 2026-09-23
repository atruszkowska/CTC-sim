#include "../../include/boundaries/round_wall.h"

/***************************************************** 
 *
 * class: Round_wall
 *
 * Class for enforcing wall boundary conditions inside
 * of a circular encolosure
 *
 * 
 ******************************************************/

/// Function for imposing BC on this boundary
void Round_wall::apply_boundary_condition(double& x, double& y)
{
	// Translate the cartesian coordinates into a circle with center at 0, 0
	x -= x_c;
	y -= y_c;	

	// Correct if not inside the allowable enclosure
	if ((x*x + y*y) > R_max*R_max)
	{
		const double theta = std::atan2(y, x);
		x = R_max*std::cos(theta);
		y = R_max*std::sin(theta);
	}

	// Translate the (possibly) modified cartesian coordinates back into any circle
	x += x_c;
	y += y_c;
}

/// Function for imposing particle-specific BC on this boundary
void Round_wall::apply_boundary_condition_for_particle(double& x, double& y, double rmin)
{
	// Translate the cartesian coordinates into a circle with center at 0, 0
	x -= x_c;
	y -= y_c;
	const double R_eff = R_max - rmin;
	// Correct if not inside the allowable enclosure
	if ((x*x + y*y) > R_eff*R_eff)
	{
		const double theta = std::atan2(y, x);
		x = R_eff*std::cos(theta);
		y = R_eff*std::sin(theta);
	}
	// Translate the (possibly) modified cartesian coordinates back into any circle
	x += x_c;
	y += y_c;
}

