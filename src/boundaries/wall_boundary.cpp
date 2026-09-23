#include "../../include/boundaries/wall_boundary.h"

/***************************************************** 
 *
 * class: Wall_boundary
 *
 * Class for enforcing wall boundary conditions
 *
 * This also acts in pairs - it assumes each wall has
 * another wall opposite to it (to be generalized)
 * 
 ******************************************************/

// Function for imposing BC on this boundary
void Wall_boundary::apply_boundary_condition(double& x)
{
	// Transform x into the 0 -> L domain world
	x -= dL;
	// Shift x if necessary
	x = std::min(std::max(x, R_min), L - R_min); 
	// Transform x back to pos_0 -> pos_f domain world
	x += dL;
}

// Function for imposing particle-specific BC on this boundary
void Wall_boundary::apply_boundary_condition_for_particle(double& x, double rmin)
{
	// Transform x into the 0 -> L domain world
	x -= dL;
	// Shift x if necessary
	x = std::min(std::max(x, rmin), L - rmin); 
	// Transform x back to pos_0 -> pos_f domain world
	x += dL;
}	
