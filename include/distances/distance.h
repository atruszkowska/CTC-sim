#ifndef DISTANCE_H
#define DISTANCE_H

#include "../common.h"

/***************************************************** 
 *
 * class: Distance
 *
 * Abstract base class for distance measurement
 * 
 ******************************************************/

class Distance {
public:
	/// Set default behavior 
	Distance() = default;
	
	/// Function for measuring distance between particles 
	/// The arguments are coordinates of particles (e.g. x1, x2 or y1, y2)
	virtual double compute_distance(const double, const double) = 0;

	/// Destructor
	virtual ~Distance() = default;
};

#endif
