#ifndef PARTICLE_H
#define PARTICLE_H

#include "../common.h"

/***************************************************** 
 *
 * class: Particle 
 * 
 ******************************************************/

class Particle {
public:

	//
	// Constructors
	//

	Particle() = default;

	Particle(const double x0, const double y0, const double theta0, const double vmag) : 
					x(x0), y(y0), theta(theta0), old_theta(theta0), new_theta(theta0), vel_mag(vmag) { compute_velocity(); needs_position_update = false; }

	Particle(const double x0, const double y0, const double theta0, const double vmag, const double r_min) : 
					x(x0), y(y0), theta(theta0), old_theta(theta0), new_theta(theta0), vel_mag(vmag), R_min(r_min) { compute_velocity(); needs_position_update = false; }

	//
	// Major functions
	//

	/// Compute particle velocity in x and y
	inline void compute_velocity() 
		{ vel_x = vel_mag*std::cos(theta); vel_y = vel_mag*std::sin(theta); }
	/// Compute new particle position
	inline void compute_new_position() 
	{
/*		if (wall_bond)
		{
			;
		} else { 
*/			x += (dx_vel + dx_bond + dx_wall_bond + dx_lift);  
			y += (dy_vel + dy_bond + dy_wall_bond + dy_lift);
//		} 
	}

	/// Contributions from external flow and/or particle motion
	inline void store_flow_contributions(const double dt)
	{ 
		dx_vel = vel_x*dt;
		dy_vel = vel_y*dt;
	}

	/// Store the newly computed heading without updating
	inline void set_new_heading(double theta_neighbors) { new_theta = theta_neighbors; }
	/// Store the newly computed position without updating
	inline void store_new_position(const std::vector<double>&& new_pos) 
	{ 
		new_x = new_pos.at(0); 
		new_y = new_pos.at(1); 
		needs_position_update = true; 
	}
	/// Set the current heading to the newly computed one
	inline void update_heading() 
	{ 
		// Update for the next angular velocity computation
		old_theta = theta;
		// Update with bond contribution
		theta = new_theta + delta_theta + delta_theta_wall_bond;
	}	
	/// Set the current position to the newly computed one
	inline void update_position() 
	{ 
		if (needs_position_update == true) 
		{
			needs_position_update = false;
			x = new_x; 
			y = new_y; 
		}
	}
	/// Compute angular velocity of the particle
	inline void compute_angular_velocity(const double delta_t) { omega = (theta - old_theta)/delta_t; }
	/// Perturb the heading with noise
	inline void add_noise(const double noise) { theta += noise; }
	/// Include external flow
	inline void add_external_flow_contribution(const double ux, const double uy) { vel_x += ux; vel_y += uy; } 
	/// Set the new position
	inline void set_new_position(const std::vector<double>&& new_pos) { x = new_pos.at(0); y = new_pos.at(1); }
	/// Add a correction to the position 
	inline void correct_x_position(const double delta_x) { x += delta_x; }
	inline void correct_y_position(const double delta_y) { y += delta_y; }
	/// Add a correction to the heading
	inline void correct_heading(const double delta_theta) { theta += delta_theta; }
	// The above, but store only for bonds, wall bonds, and lift forces (bonds are incremental)
	inline void store_bond_contribution_x(const double delta_x) { dx_bond += delta_x; }
	inline void store_bond_contribution_y(const double delta_y) { dy_bond += delta_y; }
	inline void store_bond_contribution_heading(const double _delta_theta) { delta_theta += _delta_theta; }

	inline void store_wall_bond_contribution_x(const double delta_x) { dx_wall_bond += delta_x; }
	inline void store_wall_bond_contribution_y(const double delta_y) { dy_wall_bond += delta_y; }
	inline void store_wall_bond_contribution_heading(const double _delta_theta) { delta_theta_wall_bond += _delta_theta; }

	inline void store_lift_contribution_x(const double delta_x) { dx_lift = delta_x; }
	inline void store_lift_contribution_y(const double delta_y) { dy_lift = delta_y; }

	// Increase number of bonds formed by this particle
	inline void increase_bond_number() { ++bond_number; }
	// Decrease number of bonds formed by this particle
	inline void decrease_bond_number() { --bond_number; }

	/// Particle type
	/// Two Type Identity (A=0, B=1)
	inline bool is_type(size_t t) const noexcept { return type_id_ == t; }
	inline void set_type_id(size_t t) noexcept { type_id_ = t; }
	inline size_t get_type_id() const noexcept { return type_id_; }
	inline double get_wall_clearance() const noexcept { return R_min; }	

	//
	// Getters
	//

	inline std::vector<double> get_position() const { return {x, y}; } 
	inline std::vector<double> get_all_x_contributions() const { return {dx_vel, dx_bond, dx_wall_bond, dx_lift}; }
	inline std::vector<double> get_all_y_contributions() const { return {dy_vel, dy_bond, dy_wall_bond, dy_lift}; }
	inline double get_new_heading() { return new_theta; }
	inline std::vector<double> get_heading_increments() { return {delta_theta, delta_theta_wall_bond}; }
	// These are used to set the position as well
	inline double& get_x_position() { return x; }
	inline double& get_y_position() { return y; }
	inline std::vector<double> get_velocity() const { return {vel_x, vel_y}; }
	inline double get_x_velocity() const { return vel_x; }
	inline double get_y_velocity() const { return vel_y; }
	double get_angular_velocity() const { return omega; }
	inline double get_velocity_magnitude() const { return std::sqrt(vel_x*vel_x + vel_y*vel_y);}
	inline double get_heading() const { return theta; }
	inline double get_mass() const { return mass; }
	inline double get_bond_number() const { return bond_number; }
	inline double get_radius() const { return Rp; }

	//
	// Setters
	//

 	inline void set_angular_velocity(const double _omega) { omega = _omega; } 
	inline void set_translational_velocity(const std::vector<double>& velocity) 
		{ vel_x = velocity.at(0); vel_y = velocity.at(1); }
	inline void set_particle_mass(const double _mass) { mass = _mass; } 
	inline void set_heading(const double _theta) { theta = _theta; }
	inline void set_radius(const double _Rp) { Rp = _Rp; }
	inline void reset_bond_contributions() 
	{ 
		delta_theta = 0.0; dx_bond = 0.0; dy_bond = 0.0; 
		delta_theta_wall_bond = 0.0; dx_wall_bond = 0.0; dy_wall_bond = 0.0;
	}
	inline void set_wall_bond(const bool wb) { wall_bond = wb; }
	
	// Destructor
	virtual ~Particle() = default;

private:

	//
	// Properties
	//
	
	// Position 
	double x = 0.0, y = 0.0, new_x = 0.0, new_y = 0.0;
	bool needs_position_update = false;
	double dx_vel = 0.0, dy_vel = 0.0, dx_bond = 0.0, dy_bond = 0.0; 
	double dx_wall_bond = 0.0, dy_wall_bond = 0.0; 
	double dx_lift = 0.0, dy_lift = 0.0;
	// Velocity
	double vel_x = 0.0, vel_y = 0.0;
	// Angular velocity
	double omega = 0.0;
	// Orientation
	double theta = 0.0;
	double old_theta = 0.0;
	double new_theta = 0.0;
	double delta_theta = 0.0, delta_theta_wall_bond = 0.0;
	// Speed
	double vel_mag = 0.0;
	// Mass
	double mass = 0.0;
	// Radius
	double Rp = 0.0;
	// Particle type ID
	size_t type_id_ = 0;
	// Bond count
	size_t bond_number = 0;
	// Minimum distance from the wall
	double R_min = 0.0;
	// Wall bond
	bool wall_bond = false;	

};

#endif
