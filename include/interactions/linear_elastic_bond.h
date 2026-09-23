#ifndef LINEAR_ELASTIC_BOND_H
#define LINEAR_ELASTIC_BOND_H

/***************************************************** 
 *	
 * class: Linear_elastic_bond
 * 
 * Fn is linear elastic bond, Fs and Mb is 
 * DEM-parallel bond style interaction
 * 
 ******************************************************/

#include "../common.h"
#include "../rng.h"
#include "../geometry/geometry.h"
#include "../particles/particle.h"
#include "bond.h"

class Linear_elastic_bond : public Bond {

public:

	//
	// Constructors
	//

	/// Default options 
	Linear_elastic_bond() = default;

	/// Custom options (normal and shear stiffnesses and strengths; radius; particle IDs) 
	Linear_elastic_bond(const double _k_n, const double _k_s, const double _sigma_max, 
			const double _tau_max, const double _R_bond, const size_t _pID_1, const size_t _pID_2) 
		: Bond(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, _pID_1, _pID_2) { }

	//
	// Major functions
	//

	/// Update displacements, distance is particle distance
	void add_displacement_increments(const std::unique_ptr<Particle>& particle_i, const std::unique_ptr<Particle>& particle_j, RNG& rng, const double distance, const double delta_t) override;
	/// Check for mechanical failure
	bool bond_failed(bool& failed_in_shear, bool& failed_in_tension) override;

private:
	
	// Rest separation; set on first force evaluation
    double r0 = -1.0;   
};

#endif
