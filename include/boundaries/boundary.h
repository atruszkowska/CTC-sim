#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../common.h"

/***************************************************** 
 *
 * class: Boundary
 *
 * Abstract base class for system boundaries 
 * 
 ******************************************************/

class Boundary {
public:
	/// Set default behavior 
	Boundary() = default;
	
	/// Function for imposing BC on this boundary
	/// The argument is a coordinate of a particle
	virtual void apply_boundary_condition(double&) = 0;

	/// Function for imposing BC on this boundary
	/// The arguments are coordinates of a particle
	virtual void apply_boundary_condition(double&, double&) = 0;

	/// Function for imposing particle-specific BC on this boundary
	/// The arguments are coordinate of a particle and minimum distance to the boundary
	virtual void apply_boundary_condition_for_particle(double&, double) = 0;

	/// Destructor
	virtual ~Boundary() = default;
};

#endif
