#include "../../include/boundaries/periodic_boundary.h"

/***************************************************** 
 *
 * class: Periodic_boundary
 *
 * Class for enforcing a periodic boundary condition
 * 
 ******************************************************/

// Function for imposing BC on this boundary
void Periodic_boundary::apply_boundary_condition(double& x)
{
	// Transform x into the 0 -> L domain world
	x -= dL;
	if (( x < 0 ) || ( x > L)) 
	{
		if (x < 0) {
			double temp = std::fmod((-1)*x, L);
			x = L - temp;
		} else {
			x = std::fmod(x, L);
		}
	}
	// Transform x back to pos_0 -> pos_f domain world
	x += dL;
}
