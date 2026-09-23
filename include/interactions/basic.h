#ifndef BASIC_H
#define BASIC_H

/***************************************************** 
 *
 * class: Basic
 * 
 * Simple interactions that mostly rely on removing 
 *	the overlaps between cells of the same or 
 *	different types
 *
 * The detailed description can be found in:
 * 
 ******************************************************/

#include "../common.h"
#include "../rng.h"
#include "../geometry/geometry.h"
#include "../particles/particle.h"

class Basic {

public:

	//
	// Constructors
	//

	/// Default options 
	Basic() = default;

	//
	// Major functions
	//

	/// Shifts particles that are to close to their equilibrium distance d_min
	void remove_overlaps(const std::unique_ptr<Particle>& particle_i, std::unique_ptr<Particle>& particle_j, RNG& rng, const double dsq, const double d_min);

};

#endif
