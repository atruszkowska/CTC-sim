#ifndef BOND_H
#define BOND_H

/***************************************************** 
 *	
 * class: Bond
 * 
 * DEM-parallel bond style interaction
 * 
 ******************************************************/

#include "../common.h"
#include "../rng.h"
#include "../geometry/geometry.h"
#include "../particles/particle.h"

class Bond {

public:

	//
	// Constructors
	//

	/// Default options 
	Bond() = default;

	/// Custom options (normal and shear stiffnesses and strengths; radius; particle IDs) 
	Bond(const double _k_n, const double _k_s, const double _sigma_max, 
			const double _tau_max, const double _R_bond, const size_t _pID_1, const size_t _pID_2) 
		: k_n(_k_n), k_s(_k_s), sigma_max(_sigma_max), tau_max(_tau_max), R_bond(_R_bond),
			pID_1(_pID_1),  pID_2(_pID_2)
		{ A_bond = 2*R_bond; I_bond = 2.0/3.0*R_bond*R_bond*R_bond; }

	//
	// Major functions
	//

	/// Update displacements, distance is particle distance
	virtual void add_displacement_increments(const std::unique_ptr<Particle>& particle_i, const std::unique_ptr<Particle>& particle_j, RNG& rng, const double distance, const double delta_t);
	/// Check for failure
	virtual bool bond_failed(bool& failed_in_shear, bool& failed_in_tension);

	/// Bond contributions to particle velocities
	void correct_particle_positions_and_headings(std::unique_ptr<Particle>& particle_i, std::unique_ptr<Particle>& particle_j, const double delta_t);
	/// Returns true if bond failed and is up for removal
	inline const bool failed_remove() const { return failed; } 
	/// Set the failed flag manually
	inline void set_failed(bool _failed) { failed = _failed; }
	/// Set the cumulative bond force/moment state manually (for nonzero initialization)
	inline void set_bond_state(const double _F_n, const std::vector<double>& _F_s, const double _M_b)
		{ F_n = _F_n; F_s = _F_s; M_b = _M_b; }

	//
	// Getters
	//
	/// Get the IDs of bonded particles
	inline const std::vector<size_t> get_particle_IDs() { return {pID_1, pID_2}; } 
	inline const double get_normal_force() const { return F_n; }
    inline const std::vector<double> get_shear_force() const { return F_s; }
    inline const double get_bending_moment() const { return M_b; }

    virtual ~Bond() = default;	

protected:
	
	//
	// Bond properties
	//

	// Normal stiffness of the bond
	double k_n = 0.0;
	// Shear stiffness of the bond
	double k_s = 0.0;
	// Normal strength of the bond 
	double sigma_max = 0.0;
	// Shear strength of the bond 
	double tau_max = 0.0;
	// Bond radius
	double R_bond = 0.0;
	// Bond area
	double A_bond = 0.0; 
	// Bond moment of inertia 
	double I_bond = 0.0;
	// Coordinates of the contact plane
	std::vector<double> xc{0.0, 0.0};
	// Normal of the contact plane
	std::vector<double> nc{0.0, 0.0};
	// Failed flag
	bool failed = false;

	//
	// Bonded particles
	//

	// For now IDs correspond to the position of
	// the particles in ctc_sim API pointer vector
	size_t pID_1 = 0, pID_2 = 0;

	//
	// Bond evolution
	//

	// Cumulative normal force
	double F_n = 0.0;
	// Cumulative shear force
	std::vector<double> F_s{0.0, 0.0};
	// Cumulative bending moment 
	double M_b = 0.0;
};

#endif
