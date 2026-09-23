#include "../../include/interactions/linear_elastic_bond.h"

/***************************************************** 
 *	
 * class: Linear_elastic_bond
 * 
 * Fn is linear elastic bond, Fs and Mb is 
 * DEM-parallel bond style interaction
 * 
 ******************************************************/

// Update displacements
void Linear_elastic_bond::add_displacement_increments(const std::unique_ptr<Particle>& particle_i, const std::unique_ptr<Particle>& particle_j, RNG& rng, const double dist, const double delta_t)
{
	// Particle coordinates
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
        double len_nc = std::sqrt(nc.at(0)*nc.at(0) + nc.at(1)*nc.at(1));
        nc.at(0) /= len_nc;
        nc.at(1) /= len_nc;
    }
	// Position of the contact plane
	xc.at(0) = xi + dist/2.0*nc.at(0); 
	xc.at(1) = yi + dist/2.0*nc.at(1);

	// Particle translational velocity
	std::vector<double> temp_vel = particle_i->get_velocity();
	const double vxi = temp_vel.at(0), vyi = temp_vel.at(1);
	temp_vel = particle_j->get_velocity();
	const double vxj = temp_vel.at(0), vyj = temp_vel.at(1);
	// Particle angular velocity
	const double omega_i = particle_i->get_angular_velocity();
	const double omega_j = particle_j->get_angular_velocity();

	// Velocity of the contact
	const double vcx = (vxj-omega_j*(xc.at(1)-yj))-(vxi-omega_i*(xc.at(1)-yi));
	const double vcy = (vyj+omega_j*(xc.at(0)-xj))-(vyi+omega_i*(xc.at(0)-xi));

	// Normal displacement increments
    // Capture rest length at bond formation (first evaluation)
    // rest length = contact separation, NOT formation distance          
    if (r0 < 0.0) {                                                      
        r0 = particle_i->get_radius() + particle_j->get_radius();
    }

    // Normal force: position-based overdamped spring (NOT velocity-
    // accumulated). F_n > 0 in tension (r > r0), consistent with the
    // tensile-failure sign convention in bond_failed().
    F_n = k_n * A_bond * (dist - r0);

	// Normal displacement increments
	const double delta_n = (vcx*nc.at(0) + vcy*nc.at(1))*delta_t;
	const double dun_x = delta_n*nc.at(0), dun_y = delta_n*nc.at(1); 
	// Shear displacement increments
	const double dus_x = vcx*delta_t - dun_x, dus_y = vcy*delta_t - dun_y;

	// Increment shear force
	F_s.at(0) -= (k_s*A_bond*dus_x); 
	F_s.at(1) -= (k_s*A_bond*dus_y);
	// Increment bending moment
	M_b -= (k_n*I_bond)*(omega_j - omega_i)*delta_t; 

}

// Check for failure
bool Linear_elastic_bond::bond_failed(bool& failed_in_shear, bool& failed_in_tension)
{
	// Default
	failed_in_shear = false;
	failed_in_tension = false;

	// Vector magnitude
	auto vec_mag = [](const double vx, const double vy) { return std::sqrt(vx*vx + vy*vy); };
	// Tensile stress (corrected because of a different force formulation)
	const double sigma = F_n/A_bond + std::fabs(M_b)*R_bond/I_bond;
	// Shear stress
	const double tau = vec_mag(F_s.at(0), F_s.at(1))/A_bond;

	// Failure check
	if (sigma >= sigma_max)
	{
		failed_in_tension = true;
		failed = true;
		return true;
	}
	if (tau >= tau_max)
	{
		failed_in_shear = true;
		failed = true;
		return true;
	}
	
	// Bond didn't fail
	return false;
}


