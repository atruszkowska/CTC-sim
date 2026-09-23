#include "../../include/interactions/basic.h"

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

// Shifts particles that are to close to their equilibrium distance d_min
void Basic::remove_overlaps(const std::unique_ptr<Particle>& particle_i, std::unique_ptr<Particle>& particle_j, RNG& rng, const double dsq, const double d_min)
{
	// Shift in opposite directions if closer than d_min 
	const double dist = std::sqrt(dsq);
	std::vector<double> nc(2);
	double len_nc = 0.0;

	std::vector<double> temp_coords = particle_i->get_position();
	const double xi = temp_coords.at(0); 			
	const double yi = temp_coords.at(1);
	temp_coords = particle_j->get_position();
	const double xj = temp_coords.at(0); 			
	const double yj = temp_coords.at(1);

	// Normal of the contact plane
	nc.at(0) = (xj - xi)/dist;
	nc.at(1) = (yj - yi)/dist;
	// Special case - almost fully overlapping particles
	if (std::isnan(nc.at(0) + nc.at(1))) {
		// If distance is almost zero, pick a random direction
		nc.at(0) = rng.get_random(0, 1); 	
		nc.at(1) = rng.get_random(0, 1);
		len_nc = std::sqrt(nc.at(0)*nc.at(0) + nc.at(1)*nc.at(1));
		nc.at(0) /= len_nc;
		nc.at(1) /=	len_nc;
	}

	// Shift particle j relative to i
	particle_j->store_new_position({xi + d_min*nc.at(0), yi + d_min*nc.at(1)});
}
