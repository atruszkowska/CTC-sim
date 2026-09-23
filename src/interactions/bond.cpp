#include "../../include/interactions/bond.h"

/***************************************************** 
 *	
 * class: Bond
 * 
 * DEM-parallel bond style interaction
 * 
 ******************************************************/

// Update displacements
void Bond::add_displacement_increments(const std::unique_ptr<Particle>& particle_i, const std::unique_ptr<Particle>& particle_j, RNG& rng, const double dist, const double delta_t)
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
	const double delta_n = (vcx*nc.at(0) + vcy*nc.at(1))*delta_t;
	const double dun_x = delta_n*nc.at(0), dun_y = delta_n*nc.at(1); 
	// Shear displacement increments
	const double dus_x = vcx*delta_t - dun_x, dus_y = vcy*delta_t - dun_y;

	// Increment normal force
	F_n += k_n*A_bond*delta_n;
	// Increment shear force
	F_s.at(0) -= (k_s*A_bond*dus_x); 
	F_s.at(1) -= (k_s*A_bond*dus_y);
	// Increment bending moment
	M_b -= (k_n*I_bond)*(omega_j - omega_i)*delta_t; 

}

// Check for failure
bool Bond::bond_failed(bool& failed_in_shear, bool& failed_in_tension)
{
	// Default
	failed_in_shear = false;
	failed_in_tension = false;

	// Vector magnitude
	auto vec_mag = [](const double vx, const double vy) { return std::sqrt(vx*vx + vy*vy); };
	// Tensile stress 
	const double sigma = -F_n/A_bond + std::fabs(M_b)*R_bond/I_bond;
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

// Bond contributions to particle positions and headings
void Bond::correct_particle_positions_and_headings(std::unique_ptr<Particle>& particle_i, 
					std::unique_ptr<Particle>& particle_j, const double delta_t)
{
	// Resultant on the bond directed from particle i to j
	const double F_tot_x = F_n*nc.at(0) + F_s.at(0);
	const double F_tot_y = F_n*nc.at(1) + F_s.at(1);

	// Particle masses
	const double mass_i = particle_i->get_mass(); 
	const double mass_j = particle_j->get_mass();

	// Adjust coordinates 
	particle_i->store_bond_contribution_x((delta_t*delta_t)*F_tot_x/mass_i);
	particle_i->store_bond_contribution_y((delta_t*delta_t)*F_tot_y/mass_i);

	particle_j->store_bond_contribution_x(-(delta_t*delta_t)*F_tot_x/mass_j);
	particle_j->store_bond_contribution_y(-(delta_t*delta_t)*F_tot_y/mass_j);

 	// Adjust headings
	particle_i->store_bond_contribution_heading((delta_t*delta_t)*M_b/(0.5*mass_i*R_bond*R_bond));
	particle_j->store_bond_contribution_heading(-(delta_t*delta_t)*M_b/(0.5*mass_j*R_bond*R_bond));

}

