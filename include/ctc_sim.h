#ifndef EASY_FLOW_H
#define EASY_FLOW_H

#include "ctc_sim_include.h"

#define EASY_PI 3.14159265359 

class CTC_sim 
{
 	// Function that returns a velocity value 
	// at a given position x and y
	typedef double(*vel)(double, double);

 	// Function that returns an external force value 
	// at given position x or y
	typedef double(*eforce)(double);

public:

	//
	// Constructors
	//

	/// Default initialize everything then adjust manually
	CTC_sim() = default;

	/// Introduce an external flow defined as functions of x and y, rest is default
	CTC_sim(vel v_x, vel v_y) : external_flow(v_x, v_y) { }

	/// Introduce an external flow and external force field defined as functions of x and y, rest is default
	CTC_sim(vel v_x, vel v_y, eforce F_x, eforce F_y) : external_flow(v_x, v_y), external_force(F_x, F_y) { }

	//
	// Setup
	//

	/** 
	 * \brief Define a rectangular domain - enclosure
	 * @param Lx - domain dimension in x (horizontal coordinate) 
	 * @param Ly - domain dimension in y (vertical coordinate)
	 * @param x_c - domain center, x (horizontal) coordinate
	 * @param y_c - domain center, y (vertical) coordinate
	 */
	inline void create_rectangular_enclosure(const double Lx, const double Ly, 
			const double x_c, const double y_c)
	{ geom.make_box(Lx, Ly, x_c, y_c); }

	/** 
	 * \brief Define a round (circular) domain - enclosure
	 * @param R - domain radius 
	 * @param x_c - domain center, x (horizontal) coordinate
	 * @param y_c - domain center, y (vertical) coordinate
	 */
	inline void create_circular_enclosure(const double R, const double x_c, 
		const double y_c)
	{ geom.make_circle(R, x_c, y_c); }

	/// Sets the periodic boundary conditions on all four walls
	/// First condition is for x and second is for y
	void set_all_boundaries_periodic();

	/// Sets the wall boundary conditions on all four walls
	/// First condition is for x and second is for y
	/// R_min is the minimum allowable distance with the wall
	void set_all_boundaries_walls(const double R_min);

	// Sets the wall boundary conditions on all four walls without particle radius
	void set_all_boundaries_walls();

	/// Sets LJ on all four walls
	/// First condition is for x and second is for y
	void set_all_boundaries_lj(const double sigma, const double epsilon, const double alpha, const double Fc, const int m, const int n, const double cell_mass);

	/// Sets different boundary conditions in x and y directions
	/// First argument sets x type and the second y type
	/// Current options: "periodic", "walls", or "lj-walls"
	void set_mixed_boundary_conditions(const std::string& x_type, const std::string& y_type, const double R_min = 0.0, const double sigma = 0.0, const double epsilon = 0.0, const double alpha = 0.0, const double Fc = 0.0, const int m = 0, const int n = 0, const double cell_mass = 0.0);
	
	/// Sets different boundary conditions in x and y directions for particles without minimum radius
	//void set_mixed_boundary_conditions(const std::string& x_type, const std::string& y_type, const double sigma = 0.0, const double epsilon = 0.0, const double alpha = 0.0, const double Fc = 0.0, const int m = 0, const int n = 0, const double cell_mass = 0.0);
	
	/// Set all boundary conditions to round enclosure
	void set_all_boundaries_round_wall(const double R_min);

	/// Set all boundary conditions to round enclosure without particle radius
	void set_all_boundaries_round_wall();

	/// Sets all distance computation to periodic boundaries
	void set_all_distances_periodic();

	// Sets all distance computation to wall boundaries
	void set_all_distances_walls();

	/// Sets different distance computation methods in x and y directions
	/// First argument sets x type and the second y type
	/// Current options: "periodic" or "walls"
	void set_mixed_distances(const std::string& x_type, const std::string& y_type);

	/// Sets the time step and noise levels
	inline void set_time_and_noise(const double dt, const double eta) 
		{ delta_t = dt; noise_mag = eta; }
	
	/// Sets the time step and noise levels (mixed particle types)
	inline void set_time_and_noise(const double dt, const std::vector<double> eta) 
		{ delta_t = dt; noise_magA = eta.at(0); noise_magB = eta.at(1); noise_mag = 0.5 * (noise_magA + noise_magB); }

	/// Change the range of angles for initial heading
	inline void set_initial_heading_interval(const double min_theta, const double max_theta) 
		{ min_angle = min_theta; max_angle = max_theta; }

	/// Lift force parameters
	inline void set_lift_force_coefficient(const double _CL_single, const double _CL_cluster)
		{ CL_single = _CL_single;  CL_cluster = _CL_cluster; }

	//
	// Initialization
	//

	/// Generates Np randomly heading and positioned particles with equal speed
	/// \detail x_lim, y_lim are position limits
	void initialize_particles_from_limits(const size_t Np, 
							const std::vector<double> x_lim,
							const std::vector<double> y_lim,
							const double vmag, double mass = 0.0);

	/// Positions particles at x y locations specified by the x/y vectors 	
	/// Heading is random
	void initialize_particles_from_list(const size_t Np, 
							const std::vector<double>& x_pos,
							const std::vector<double>& y_pos,
							const double vmag, double mass = 0.0);

	/// Generates Np randomly heading and positioned particles with equal speed and multiple different types
	/// \detail x_lim, y_lim are position limits
	void initialize_particles_from_limits(const std::vector<size_t> Np_list,
							const std::vector<double> x_lim,
							const std::vector<double> y_lim,
							const std::vector<double> vmag,
							const std::vector<size_t> type_id,
							const std::vector<double> R_min, double mass = 0.0);
	

	/// Generates Np randomly heading and positioned particles inside of a circular enclosure
	/// \detail center_coordinates - coordinates of the enclosure center
	void initialize_particles_round_enclosure(const size_t Np, 
							const std::vector<double> center_coordinates,
							double enclosure_radius, const double vmag, double mass = 0.0);

	/// Generates Np randomly heading and positioned particles inside of a circular enclosure with two different types
	/// \detail center_coordinates - coordinates of the enclosure center
	void initialize_particles_round_enclosure(const std::vector<size_t> Np,
							const std::vector<double> center_coordinates,
							double enclosure_radius, const std::vector<double> vmag,
							const std::vector<size_t> type_id, 
							const std::vector<double> R_min, double mass = 0.0);

	/// Re-initializes particles after deleting pointer 
	void reinitialize_particles(const size_t Np,
						const std::vector<double> x_lim,
						const std::vector<double> y_lim,
						const double vmag, double mass = 0.0);

	/// Re-initialize particles at x y locations specified by the x/y vectors 	
	void reinitialize_particles_from_list(const size_t Np, 
							const std::vector<double>& x_pos,
							const std::vector<double>& y_pos,
							const double vmag, double mass = 0.0);

	/// Re-initializes particles after deleting pointer with multiple different types					
	void reinitialize_particles(const std::vector<size_t> Np, 
						const std::vector<double> x_lim,
						const std::vector<double> y_lim,
						const std::vector<double> vmag, 
						const std::vector<size_t> type_id,
						const std::vector<double> R_min, 
						double mass = 0.0); 

	/// Re-initializes particles in a circular enclosure
	void reinitialize_particles(const size_t Np, 
							const std::vector<double> center_coordinates,
							double enclosure_radius, const double vmag, double mass = 0.0);

	/// Re-initializes particles in a circular enclosure with multiple different types
	void reinitialize_particles(const std::vector<size_t> Np, 
							const std::vector<double> center_coordinates,
							double enclosure_radius, const std::vector<double> vmag, 
							const std::vector<size_t> type_id,
							const std::vector<double> R_min, 
							double mass = 0.0);
	

	/// Find and initialize cell bonds
	/// @param dist_to_bon - bonds will be made between cells at that and smaller distance
	size_t install_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond, 
			const double max_bond_len = 1e3, const std::string type = "DEM");

	/// Find and initialize bonds between walls and cells
	/// @param dist_to_bon - bonds will be made between cells at that and smaller distance
	size_t install_wall_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond, 
			const double max_wall_bond_len = 1e3);

	//
	// Simulation
	//
	
	/// Invoke periodic boundary conditions
	inline void apply_periodic_boundary_conditions()
		{ apply_boundary_conditions(); }

	/// Invoke wall boundary conditions
	inline void apply_wall_boundary_conditions()
		{ apply_boundary_conditions(); }
	

	/// Invoke mixed boundary conditions
	inline void apply_mixed_boundary_conditions()
		{ apply_boundary_conditions(); }

	/// Invoke xy boundary conditions (i.e. with two arguments)
	inline void apply_round_wall_boundary_conditions()
		{ apply_xy_boundary_conditions(); }

	/// Invoke wall boundary conditions for multiple particle types
	inline void apply_wall_boundary_conditions_for_particle()
		{ apply_boundary_conditions_for_particle(); }
	/// Invoke round wall boundary conditions for multiple particle types
	inline void apply_round_wall_boundary_conditions_for_particle()
		{ apply_xy_boundary_conditions_for_particle(); }
	/// Invoke mixed boundary conditions for multiple particle types
	inline void apply_mixed_boundary_conditions_for_particle()
		{ apply_boundary_conditions_for_particle(); }


	/// Apply external force field (displaces particles)
	void apply_external_force_field();

	/// Find the average directions of the neighbors (particles within radius R)
	void compute_average_heading(const double R);	

	/// Find the average directions of the neighbors (particles within radius R) in the presence of an external flow
	void compute_average_heading(const double R, const bool needs_external_flow);

	/// Find the average directions of the neighbors (particles within radius R) for two different particle types
	void compute_average_heading_two_types(const std::vector<double> R, const bool needs_external_flow);

	/// Move particles following their current velocity
	void propagate_particles();
	
	/// Move particles following their current velocity and external flow
	void propagate_particles(const bool needs_external_flow);

	/// Move away particles closer than d_min
	void remove_particle_overlaps(const double d_min); 

	/// Move away particles closer than d_min for two different particle types
	void remove_particle_overlaps_two_types(const std::vector<double> d_min);

	/// Compute bond forces and moments, check failure, correct particle properties
	void process_bonds();

	/// Compute wall bond forces and moments, check failure, correct particle properties
	void process_wall_bonds();

	/// Compute and store lift force contribution for each particle
	void add_lift_force(const double Rp, const double rho, const double mu);

	/// Install new bonds or reinstall old one if within proximity limits
	void reinstall_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond);

	/// Install new wall bonds or reinstall old one if within proximity limits
	void reinstall_wall_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond);

	/// Save the number of bonds at this time step
	inline void track_bonds() { bond_number.push_back(particle_bonds.size()); }

	/// Initialize all cell bonds with nonzero force/moment state
	/// Applies the same (F_n0, F_s0, M_b0) to every existing bond.
	/// Call after install_wall_bonds. F_s0 is the shear force vector {x, y}.
	void initialize_bond_state(const double F_n0,
			const std::vector<double>& F_s0, const double M_b0);

	/// Initialize all wall bonds with nonzero force/moment state
	/// Applies the same (F_n0, F_s0, M_b0) to every existing wall bond.
	/// Call after install_wall_bonds. F_s0 is the shear force vector {x, y}.
	void initialize_wall_bond_state(const double F_n0,
			const std::vector<double>& F_s0, const double M_b0);
		
	/// Compute the new positions and headings from all contributions
	void update_positions_and_headings();
	
	//
	// Getters
	// 
	
	/// Returns the lower and upper limits of a rectangular enclosure [x0, y0, xf, yf]
    inline std::vector<double> get_enclosure_limits() const { return geom.get_box_limits(); }

	/// Returns round enclosure radius
    inline double get_enclosure_radius() const { return geom.get_radius(); }

	/// Returns round enclosure center coordinates [x, y]
    inline std::vector<double> get_enclosure_center() const { return geom.get_center(); }

	/// Constant reference to a vector of unique pointers to all particles
	inline const std::vector<std::unique_ptr<Particle>>& get_particles() const { return particles; }

	/// Returns a vector with bond number at a specific time
	inline std::vector<size_t> get_bond_number_with_time() const { return bond_number; }

	/// Return failure mode vector	
	// (0 - bond failed in tension, 1 - in shear, 2 - bond too long)
	inline std::vector<std::vector<int>> get_failure_mode_data() const { return bond_failure_mode; }

	inline std::vector<std::vector<int>> get_wall_failure_mode_data() const { return wall_bond_failure_mode; }

	/// IDs of currently bonded cells: lower ID-higher ID
	inline const std::set<std::string>& get_bonded_pairs_IDs() const { return bonded_pairs; }

	/// IDs of currently bonded cells: lower ID-higher ID
	inline const std::set<std::string>& get_wall_bonded_pairs_IDs() const { return wall_bonded_pairs; }

	/// Bonds
	inline std::list<std::unique_ptr<Bond>>& get_particle_bonds() { return particle_bonds; }
	inline std::list<Wall_bond>& get_wall_bonds() { return wall_bonds; }

	//
	// Postprocessing
	//

	/// Mean velocity magnitude over all particles - not normalized
	double compute_average_velocity_magnitude();
	/// Polarization
	double compute_polarization();

	// 
	// Testing utilities
	// 
	
	/// Install a bond between particle pair
	void create_a_bond(const size_t ID_1, const size_t ID_2, 
			const double _k_n, const double _k_s, const double _sigma_max,
            const double _tau_max, const double _R_bond);

	/// Remove a bond between particle pair
	void remove_a_bond(const size_t ID_1, const size_t ID_2);
	

private:

	//
	// Members
	//

	// Time step
	double delta_t = 0.0;

	// Noise level
	double noise_mag = 0.0;
	double noise_magA = 0.0;
	double noise_magB = 0.0;

	// BC flags
	std::string x_bc_type, y_bc_type;

	// Initial heading limits
	double min_angle = -EASY_PI;
	double max_angle = EASY_PI;

	// Other parameters
	double CL_single = 0.0, CL_cluster = 0.0;
	double max_bond_dist = 0.0, max_wall_bond_dist = 0.0;

	// Random distribution generator
	RNG rng;
	
	// Physical geometry of the computational domain
	Geometry geom;
	
	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// Distance computation
	// First instance is the approach to compute the distance in x
	// Second instance is the approach to compute the distance in y
	std::vector<std::unique_ptr<Distance>> distances;

	// Boundaries
	std::vector<std::unique_ptr<Boundary>> boundaries;
	std::vector<std::unique_ptr<Boundary_xy>> xy_boundaries;

	// Interactions
	Basic basic_interactions;

	// External flow field
	Flow external_flow;	

	// External force field
	Force external_force; 

	// Type of particle bonds
	std::string particle_bond_type;
	// All current bonds between particles
	std::list<std::unique_ptr<Bond>> particle_bonds;
	// Bonded particle pairs (lower particle ID-higher particle ID)
	std::set<std::string> bonded_pairs;
	// Vector with bond numbers at user-specified times
	std::vector<size_t> bond_number;
	// Vector with number of failed bonds numbers
	std::vector<size_t> failed_bond_numbers;
	// Failure mode vector 
	// (0 - bond failed in tension, 1 - in shear, 2 - bond too long)
	// Each inner vector are bonds that failed at that particular time
	std::vector<std::vector<int>> bond_failure_mode; 

	// All current bonds between particles and walls
	std::list<Wall_bond> wall_bonds;
	// Bonded particle-wall pairs (lower particle ID-wall ID placeholder)
	std::set<std::string> wall_bonded_pairs;
	// Vector with bond numbers at user-specified times
	std::vector<size_t> wall_bond_number;
	// Failure mode vector (0 - bond failed in tension, 1 - in shear)
	// Each inner vector are bonds that failed at that particular time
	std::vector<std::vector<int>> wall_bond_failure_mode;

	//
	// Methods
	//

	// Measure the distance between two x coordinates 
	inline double compute_x_distance(const double xi, const double xj)
		{ return (distances.at(0))->compute_distance(xi, xj); }
	// Measure the distance between two y coordinates 
	inline double compute_y_distance(const double yi, const double yj)
		{ return (distances.at(1))->compute_distance(yi, yj); }

	// Invoke boundary conditions
	void apply_boundary_conditions();
	// Invoke xy boundary conditions (i.e. with two arguments)
	void apply_xy_boundary_conditions();
	// Invoke boundary conditions for multiple particle types
	void apply_boundary_conditions_for_particle();
	// Invoke xy boundary conditions for multiple particle types
	void apply_xy_boundary_conditions_for_particle();
};

#endif
