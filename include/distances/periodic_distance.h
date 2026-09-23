#ifndef PERIODIC_DISTANCE_H
#define PERIODIC_DISTANCE_H

#include "distance.h"
#include "../common.h"

/***************************************************** 
 *
 * class: Periodic_distance
 *
 * Class for computing the distance between particles
 * across periodic boundaries.  
 * 
 ******************************************************/

class Periodic_distance : public Distance 
{
public:

	/// Sets variables needed for the computation
 	/// @param _L - distance between pairs of boundaries
	Periodic_distance(const double _L) : L(_L) { }
	
	/// Function for measuring distance between particles 
	/// The arguments are coordinates of particles (e.g. x1, x2 or y1, y2)
	inline double compute_distance(const double xi, const double xj) override
	{ return std::min<double>(std::fabs(xi - xj), L - std::fabs(xi - xj)); }

protected:

	// Length in the direction periodic boundaries are imposed
	double L = 0.0;
};

#endif
