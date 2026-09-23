#ifndef BOUNDARY_XY_H
#define BOUNDARY_XY_H

#include "../common.h"

/***************************************************** 
 *
 * class: Boundary_xy
 *
 * Abstract base class for boundaries that take two
 * arguments (coordinates, likely)
 * 
 ******************************************************/

class Boundary_xy 
{

public:

	/// Set default behavior 
	Boundary_xy() = default;
	
	/// Function for imposing BC on this boundary
	/// The arguments are x and y coordinates of a particle
	virtual void apply_boundary_condition(double&, double&) = 0;

	/// Function for imposing particle-specific BC on this boundary
	/// The arguments are coordinates of a particle and minimum distance to the boundary
	virtual void apply_boundary_condition_for_particle(double&, double&, double) = 0;
	
	/// Destructor
	virtual ~Boundary_xy() = default;

};

#endif
