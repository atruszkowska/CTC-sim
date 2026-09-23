#ifndef WALL_BOND_H
#define WALL_BOND_H

/***************************************************** 
 *	
 * class: Wall_bond
 * 
 * Fn is linear elastic bond, Fs and Mb is 
 * DEM-parallel bond style interaction 
 * 
 ******************************************************/

#include "../common.h"
#include "../rng.h"
#include "../geometry/geometry.h"
#include "../particles/particle.h"

class Wall_bond {

public:

	//
	// Constructors
	//

	/// Default options 
	Wall_bond() = default;

	/// Custom options (normal and shear stiffnesses and strengths; radius; particle ID, wall normal vector) 
	Wall_bond(const double _k_n, const double _k_s, const double _sigma_max, 
			const double _tau_max, const double _R_bond, const size_t _pID_1, const std::vector<double> wall_xy) 
		: k_n(_k_n), k_s(_k_s), sigma_max(_sigma_max), tau_max(_tau_max), R_bond(_R_bond),
			pID_1(_pID_1),  wall_x(wall_xy.at(0)), wall_y(wall_xy.at(1))
		{ A_bond = 2*R_bond; I_bond = 2.0/3.0*R_bond*R_bond*R_bond; }

	//
	// Major functions
	//

	/// Get the ID of bonded particle
	const size_t get_particle_ID() const { return pID_1; } 
	/// Update displacements, distance is particle distance
	void add_displacement_increments(const std::unique_ptr<Particle>& particle_i, RNG& rng, const double distance, const double delta_t);
	/// Check for failure
	bool bond_failed(bool& failed_in_shear, bool& failed_in_tension);
	/// Bond contributions to particle velocities
	void correct_particle_positions_and_headings(std::unique_ptr<Particle>& particle_i, const double delta_t);
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

	inline const double get_normal_force() const { return F_n; }
    inline const std::vector<double> get_shear_force() const { return F_s; }
    inline const double get_bending_moment() const { return M_b; }
	inline const double get_x_coord() const { return wall_x; }
	inline const double get_y_coord() const { return wall_y; }

private:
	
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
	// For now ID correspond to the position of
	// the particle in ctc_sim API pointer vector
	size_t pID_1 = 0;
	// Wall x coordinate (original)
	double wall_x = 0.0;
	// Wall y coordinate (original)
	double wall_y = 0.0;
	// Rest separation; set on first force evaluation
    double r0 = -1.0;   

	// Coordinates of the contact plane
	std::vector<double> xc{0.0, 0.0};
	// Normal of the contact plane
	std::vector<double> nc{0.0, 0.0};
	// Failed flag
	bool failed = false;

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
