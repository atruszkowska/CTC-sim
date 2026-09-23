#ifndef CARTESIAN_DISTANCE_H
#define CARTESIAN_DISTANCE_H

#include "distance.h"
#include "../common.h"

/***************************************************** 
 *
 * class: Cartesian_distance
 *
 * Class for computing the distance between particles
 * across wall boundaries.  
 * 
 ******************************************************/

class Cartesian_distance : public Distance 
{
public:

	/// Default behavior - no need for extra parameters 
	Cartesian_distance() = default;
	
	/// Function for measuring distance between particles 
	/// The arguments are coordinates of particles (e.g. x1, x2 or y1, y2)
	inline double compute_distance(const double x1, const double x2) override
		{ return (x1 - x2); }

};

#endif
