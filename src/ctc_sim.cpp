#include "../include/ctc_sim.h"

// Sets the periodic boundary conditions on all four walls
void CTC_sim::set_all_boundaries_periodic()
{	
	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();
	std::vector<double> all_lim = geom.get_box_limits(); 

	boundaries.push_back(std::unique_ptr<Periodic_boundary>(new Periodic_boundary(Lx, {all_lim.at(0), all_lim.at(2)})));	
	boundaries.push_back(std::unique_ptr<Periodic_boundary>(new Periodic_boundary(Ly, {all_lim.at(1), all_lim.at(3)})));
}

// Sets the wall boundary conditions on all four walls
void CTC_sim::set_all_boundaries_walls(const double R_min)
{	
	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();
	std::vector<double> all_lim = geom.get_box_limits(); 

	boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Lx, {all_lim.at(0), all_lim.at(2)}, R_min)));	
	boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Ly, {all_lim.at(1), all_lim.at(3)}, R_min)));
}

// Sets the wall boundary conditions on all four walls without particle radius
void CTC_sim::set_all_boundaries_walls()
{	
	//boundaries.clear();
	x_bc_type = "walls";
        y_bc_type = "walls";
	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();
	std::vector<double> all_lim = geom.get_box_limits(); 

	boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Lx, {all_lim.at(0), all_lim.at(2)})));	
	boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Ly, {all_lim.at(1), all_lim.at(3)})));
}

// Sets the LJ boundary conditions on all four walls
void CTC_sim::set_all_boundaries_lj(const double sigma, const double epsilon, const double alpha, const double Fc, const int m, const int n, const double cell_mass)
{
	boundaries.clear();

	x_bc_type = "lj-walls";
	y_bc_type = "lj-walls";

	const double Lx = geom.get_x_length();
    	const double Ly = geom.get_y_length();
    	std::vector<double> all_lim = geom.get_box_limits();

    	boundaries.push_back(std::unique_ptr<Lj_wall>(new Lj_wall(Lx, {all_lim.at(0), all_lim.at(2)}, sigma, epsilon, alpha, Fc, m, n, delta_t, cell_mass)));
    	boundaries.push_back(std::unique_ptr<Lj_wall>(new Lj_wall(Ly, {all_lim.at(1), all_lim.at(3)}, sigma, epsilon, alpha, Fc, m, n, delta_t, cell_mass)));

}

// Sets different boundary conditions in x and y directions
void CTC_sim::set_mixed_boundary_conditions(const std::string& x_type, const std::string& y_type, const double R_min, const double sigma, const double epsilon, const double alpha, const double Fc, const int m, const int n, const double cell_mass)
{
	boundaries.clear();

	x_bc_type = x_type;
	y_bc_type = y_type;

	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();
	std::vector<double> all_lim = geom.get_box_limits();

	if (x_type == "periodic")
	{
		boundaries.push_back(std::unique_ptr<Periodic_boundary>(new Periodic_boundary(Lx, {all_lim.at(0), all_lim.at(2)})));
	} else if (x_type == "walls") {
		boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Lx, {all_lim.at(0), all_lim.at(2)}, R_min)));
	} else if (x_type == "lj-walls") {
	    boundaries.push_back(std::unique_ptr<Lj_wall>(new Lj_wall(Lx, {all_lim.at(0), all_lim.at(2)}, sigma, epsilon, alpha, Fc, m, n, delta_t, cell_mass)));
	} else {
		throw std::invalid_argument("Wrong type of x-boundary: " + x_type);
	}

	if (y_type == "periodic")
	{
		boundaries.push_back(std::unique_ptr<Periodic_boundary>(new Periodic_boundary(Ly, {all_lim.at(1), all_lim.at(3)})));
	} else if (y_type == "walls") {
		boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Ly, {all_lim.at(1), all_lim.at(3)}, R_min)));
	} else if (y_type == "lj-walls") {
        boundaries.push_back(std::unique_ptr<Lj_wall>(new Lj_wall(Ly, {all_lim.at(1), all_lim.at(3)}, sigma, epsilon, alpha, Fc, m, n, delta_t, cell_mass)));
	} else {
		throw std::invalid_argument("Wrong type of y-boundary: " + y_type);
	}
	
}

// Sets different boundary conditions in x and y directions for particles without minimum radius
/*void CTC_sim::set_mixed_boundary_conditions(const std::string& x_type, const std::string& y_type, const double sigma, const double epsilon, const double alpha, const double Fc, const int m, const int n, const double cell_mass)
{
	boundaries.clear();

	x_bc_type = x_type;
	y_bc_type = y_type;

	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();
	std::vector<double> all_lim = geom.get_box_limits();

	if (x_type == "periodic")
	{
		boundaries.push_back(std::unique_ptr<Periodic_boundary>(new Periodic_boundary(Lx, {all_lim.at(0), all_lim.at(2)})));
	} else if (x_type == "walls") {
		boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Lx, {all_lim.at(0), all_lim.at(2)})));
	} else if (x_type == "lj-walls") {
	    boundaries.push_back(std::unique_ptr<Lj_wall>(new Lj_wall(Lx, {all_lim.at(0), all_lim.at(2)}, sigma, epsilon, alpha, Fc, m, n, delta_t, cell_mass)));
	} else {
		throw std::invalid_argument("Wrong type of x-boundary: " + x_type);
	}
	if (y_type == "periodic")
	{
		boundaries.push_back(std::unique_ptr<Periodic_boundary>(new Periodic_boundary(Ly, {all_lim.at(1), all_lim.at(3)})));
	} else if (y_type == "walls") {
		boundaries.push_back(std::unique_ptr<Wall_boundary>(new Wall_boundary(Ly, {all_lim.at(1), all_lim.at(3)})));
	} else if (y_type == "lj-walls") {
		boundaries.push_back(std::unique_ptr<Lj_wall>(new Lj_wall(Ly, {all_lim.at(1), all_lim.at(3)}, sigma, epsilon, alpha, Fc, m, n, delta_t, cell_mass)));
	} else {
		throw std::invalid_argument("Wrong type of y-boundary: " + y_type);
	}	
	
}
*/

// Set all boundary conditions to round enclosure
void CTC_sim::set_all_boundaries_round_wall(const double R_min)
{
	xy_boundaries.push_back(std::unique_ptr<Round_wall>(new Round_wall(geom.get_radius(), geom.get_center(), R_min)));
}

// Set all boundary conditions to round enclosure without particle radius
void CTC_sim::set_all_boundaries_round_wall()
{
	xy_boundaries.push_back(std::unique_ptr<Round_wall>(new Round_wall(geom.get_radius(), geom.get_center())));
}

// Sets all distance computation to periodic boundaries
void CTC_sim::set_all_distances_periodic()
{
	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();

	distances.push_back(std::unique_ptr<Periodic_distance>(new Periodic_distance(Lx)));	
	distances.push_back(std::unique_ptr<Periodic_distance>(new Periodic_distance(Ly)));
}

// Sets all distance computation to wall boundaries
void CTC_sim::set_all_distances_walls()
{
	distances.push_back(std::unique_ptr<Cartesian_distance>(new Cartesian_distance));	
	distances.push_back(std::unique_ptr<Cartesian_distance>(new Cartesian_distance));
}

// Sets different distance computation methods in x and y directions
void CTC_sim::set_mixed_distances(const std::string& x_type, const std::string& y_type)
{
	const double Lx = geom.get_x_length();	
	const double Ly = geom.get_y_length();

	if (x_type == "periodic")
	{
		distances.push_back(std::unique_ptr<Periodic_distance>(new Periodic_distance(Lx)));	
	} else if (x_type == "walls") {
		distances.push_back(std::unique_ptr<Cartesian_distance>(new Cartesian_distance));
	} else if (x_type == "lj-walls") {
		distances.push_back(std::unique_ptr<Cartesian_distance>(new Cartesian_distance));
	} else {
		throw std::invalid_argument("Wrong type of x-distance: " + x_type);
	}

	if (y_type == "periodic")
	{
		distances.push_back(std::unique_ptr<Periodic_distance>(new Periodic_distance(Ly)));
	} else if (y_type == "walls") {
		distances.push_back(std::unique_ptr<Cartesian_distance>(new Cartesian_distance));
	} else if (y_type == "lj-walls") {
		distances.push_back(std::unique_ptr<Cartesian_distance>(new Cartesian_distance));
	} else {
		throw std::invalid_argument("Wrong type of y-distance: " + y_type);
	}
}

// Generates Np randomly heading and positioned particles
void CTC_sim::initialize_particles_from_limits(const size_t Np, 
							const std::vector<double> x_lim,
							const std::vector<double> y_lim,
							const double vmag, double mass)
{		
	double x0 = 0.0, y0 = 0.0;
	double theta = 0.0;
	for (size_t i = 0; i < Np; ++i) 
	{
		x0 = rng.get_random(x_lim.at(0), x_lim.at(1));
		y0 = rng.get_random(y_lim.at(0), y_lim.at(1));
		theta = rng.get_random(min_angle, max_angle);
		particles.push_back(std::unique_ptr<Particle>(new Particle(x0, y0, theta, vmag)));	
	}
	// Set particle masses (default is 0.0) 
	std::for_each(particles.begin(), particles.end(), [&mass](std::unique_ptr<Particle>& particle){ particle->set_particle_mass(mass); });
}

// Positions particles at x y locations specified by the x/y vectors 	
// Heading is random
void CTC_sim::initialize_particles_from_list(const size_t Np, 
							const std::vector<double>& x_pos,
							const std::vector<double>& y_pos,
							const double vmag, double mass)
{		
	double x0 = 0.0, y0 = 0.0;
	double theta = 0.0;
	for (size_t i = 0; i < Np; ++i) 
	{
		x0 = x_pos.at(i);
		y0 = y_pos.at(i);
		theta = rng.get_random(min_angle, max_angle);
		particles.push_back(std::unique_ptr<Particle>(new Particle(x0, y0, theta, vmag)));	
	}
	// Set particle masses (default is 0.0) 
	std::for_each(particles.begin(), particles.end(), [&mass](std::unique_ptr<Particle>& particle){ particle->set_particle_mass(mass); });
}

// Generates Np randomly heading and positioned particles (mixed type)
void CTC_sim::initialize_particles_from_limits(const std::vector<size_t> Np_types, const std::vector<double> x_lim,
												const std::vector<double> y_lim,
												const std::vector<double> vmag, 
												const std::vector<size_t> type_id,
												const std::vector<double> R_min, double mass)
{
	double x0 = 0.0, y0 = 0.0;
    	double theta = 0.0;
	for (size_t i = 0; i < type_id.size(); ++i)
	{
		for (size_t j = 0; j < Np_types.at(i); ++j)
		{
			x0 = rng.get_random(x_lim.at(0), x_lim.at(1));
			y0 = rng.get_random(y_lim.at(0), y_lim.at(1));
			theta = rng.get_random(min_angle, max_angle);
			auto p = std::unique_ptr<Particle>(new Particle(x0, y0, theta, vmag.at(i), R_min.at(i)));

			// Set type label (A=0, B=1)
			p->set_type_id(type_id.at(i));
			particles.push_back(std::move(p));
				
		}
	}
	// Set particle masses (default is 0.0)
	std::for_each(particles.begin(), particles.end(), [&mass](std::unique_ptr<Particle>& particle){ particle->set_particle_mass(mass); });
}

// Generates Np randomly heading and positioned particles inside of a circular enclosure
void CTC_sim::initialize_particles_round_enclosure(const size_t Np, 
							const std::vector<double> center_coordinates,
							double enclosure_radius, const double vmag, double mass)
{		
	double R0 = 0.0, x0 = 0.0, y0 = 0.0;
	double theta = 0.0;
	for (size_t i = 0; i < Np; ++i) 
	{
		R0 = rng.get_random(0.0, enclosure_radius);
		theta = rng.get_random(min_angle, max_angle);
		x0 = center_coordinates.at(0) + R0*std::cos(theta);
		y0 = center_coordinates.at(1) + R0*std::sin(theta);
		particles.push_back(std::unique_ptr<Particle>(new Particle(x0, y0, theta, vmag)));	
	}
	// Set particle masses (default is 0.0) 
	std::for_each(particles.begin(), particles.end(), [&mass](std::unique_ptr<Particle>& particle){ particle->set_particle_mass(mass); });
}

// Generates Np randomly heading and positioned particles inside of a circular enclosure (mixed type)
void CTC_sim::initialize_particles_round_enclosure(const std::vector<size_t> Np_types,
							const std::vector<double> center_coordinates,
							double enclosure_radius, const std::vector<double> vmag, 
							const std::vector<size_t> type_id, const std::vector<double> R_min, double mass)
{
		double R0 = 0.0, x0 = 0.0, y0 = 0.0;
		double theta = 0.0;
		for (size_t i = 0; i < type_id.size(); ++i) {
			for (size_t j = 0; j < Np_types.at(i); ++j)
			{
					R0 = rng.get_random(0.0, enclosure_radius);
					theta = rng.get_random(min_angle, max_angle);
					x0 = center_coordinates.at(0) + R0*std::cos(theta);
					y0 = center_coordinates.at(1) + R0*std::sin(theta);
					auto p = std::unique_ptr<Particle>(new Particle(x0, y0, theta, vmag.at(i), R_min.at(i)));

					// Set type label (A=0, B=1, ...)
					p->set_type_id(type_id.at(i));
					particles.push_back(std::move(p));
			}
				
		}
		// Set particle masses (default is 0.0)
		std::for_each(particles.begin(), particles.end(), [&mass](std::unique_ptr<Particle>& particle){ particle->set_particle_mass(mass); });


}

// Re-initializes particles 
void CTC_sim::reinitialize_particles(const size_t Np,
						const std::vector<double> x_lim,
						const std::vector<double> y_lim,
						const double vmag, double mass) 
{
	particles.clear();
	particles.reserve(Np);
	initialize_particles_from_limits(Np, x_lim, y_lim, vmag, mass);
}

// Re-initializes particles 
void CTC_sim::reinitialize_particles_from_list(const size_t Np,
							const std::vector<double>& x_pos,
							const std::vector<double>& y_pos,
							const double vmag, double mass) 
{
	particles.clear();
	particles.reserve(Np);
	initialize_particles_from_list(Np, x_pos, y_pos, vmag, mass);
}

// Re-initializes particles (mixed type)
void CTC_sim::reinitialize_particles(const std::vector<size_t> Np, 
						const std::vector<double> x_lim,
						const std::vector<double> y_lim,
						const std::vector<double> vmag, 
						const std::vector<size_t> type_id,
						const std::vector<double> R_min, 
						double mass) 
{
	particles.clear();
	particles.reserve(std::accumulate(Np.begin(), Np.end(), 0));
	initialize_particles_from_limits(Np, x_lim, y_lim, vmag, type_id, R_min, mass);	
				
}

// Re-initializes particles in a circular enclosure
void CTC_sim::reinitialize_particles(const size_t Np, 
							const std::vector<double> center_coordinates,
							double enclosure_radius, const double vmag, double mass) 
{
	particles.clear();
	particles.reserve(Np);
	initialize_particles_round_enclosure(Np, center_coordinates, enclosure_radius, vmag, mass);
}

// Re-initializes particles in a circular enclosure (mixed type)
void CTC_sim::reinitialize_particles(const std::vector<size_t> Np, 
							const std::vector<double> center_coordinates,
							double enclosure_radius, 
							const std::vector<double> vmag, 
							const std::vector<size_t> type_id,
							const std::vector<double> R_min, 
							double mass) 
{
	particles.clear();
	particles.reserve(std::accumulate(Np.begin(), Np.end(), 0));
	initialize_particles_round_enclosure(Np, center_coordinates, enclosure_radius, vmag, type_id, R_min, mass);
}

// Find and initialize cell bonds
size_t CTC_sim::install_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond, 
			const double max_bond_len, const std::string type)
{
	// Global bond type
	particle_bond_type = type;		

	// Empty all existing bond-related structures
    particle_bonds.clear();
    bond_number.clear();
    bond_failure_mode.clear();
    bonded_pairs.clear();
	
 	// Temporary variables
	std::vector<double> temp_coords;
	std::unique_ptr<Bond> temp_bond;
	Particle temp_particle_i, temp_particle_j;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;

	// Store max allowable bond length
	max_bond_dist = max_bond_len;	

	// Compute the distance between particles
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		for (size_t jp = 0; jp < particles.size(); ++jp) 
		{
		
			// Skip itself or if this pair was already considered	
			if (ip >= jp) 
			{
				continue;
			}

			// Compute interparticle distance
			temp_particle_i = *(particles.at(ip));
			temp_particle_j = *(particles.at(jp));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);
	
			dsq = dx*dx + dy*dy;

			// If particles are close enough, form a bond 
			if (dsq <= dist_to_bond*dist_to_bond) 
			{
				std::string pair = std::to_string(ip) + "-" + std::to_string(jp);
				bonded_pairs.insert(pair);
				
				if (particle_bond_type == "DEM") {
				    temp_bond = std::make_unique<Bond>(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ip, jp);
				} 
				else if (particle_bond_type == "linear") {
				    temp_bond = std::make_unique<Linear_elastic_bond>(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ip, jp);
				} 
				else {
				    throw std::invalid_argument("Wrong type of particle bond: " + particle_bond_type);
				}
				particle_bonds.push_back(std::move(temp_bond));

				// Increase the number of particle bonds
				particles.at(ip)->increase_bond_number();
				particles.at(jp)->increase_bond_number();
			}
 		}
	}
	return particle_bonds.size();
}

// Find and initialize cell-wall bonds
size_t CTC_sim::install_wall_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond,  
			const double max_wall_bond_len)
{
	// Empty all existing wall bond-related structures
    wall_bonds.clear();
    wall_bond_number.clear();
    wall_bond_failure_mode.clear();
    wall_bonded_pairs.clear();
	
 	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;
	// Coordinates of walls
	std::vector<double> walls = geom.get_box_limits();
	std::vector<double> y_coords{1,3};

	// Store max allowable bond length
	max_wall_bond_dist = max_wall_bond_len;

	// Compute the distance between particles and the wall
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		for (size_t jp = 0; jp < 2; ++jp) 
		{
			// Compute the distance
			temp_particle_i = *(particles.at(ip));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);

			xj = temp_coords.at(0); 			
			yj = walls.at(y_coords.at(jp)); 
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);
	
			dsq = dx*dx + dy*dy;

			// If particles are close enough, form a bond 
			if (dsq <= dist_to_bond*dist_to_bond) 
			{
				std::string pair = std::to_string(ip) + "-wall";
				wall_bonded_pairs.insert(pair);
			
				Wall_bond temp_bond(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ip, {xj, yj});
				wall_bonds.push_back(temp_bond);

				particles.at(ip)->set_wall_bond(true); 
			}
 		}
	}
	return wall_bonds.size();
}

// Install new bonds or reinstall old one if within proximity limits
void CTC_sim::reinstall_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond)
{

 	// Temporary variables
	std::vector<double> temp_coords;
	std::unique_ptr<Bond> temp_bond;
	Particle temp_particle_i, temp_particle_j;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;

	// Compute the distance between particles
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		for (size_t jp = 0; jp < particles.size(); ++jp) 
		{
		
			// Skip itself or if this pair was already considered	
			if (ip >= jp) 
			{
				continue;
			}

			// Compute interparticle distance
			temp_particle_i = *(particles.at(ip));
			temp_particle_j = *(particles.at(jp));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);
	
			dsq = dx*dx + dy*dy;

			// If particles are close enough and the bond is not present, form a bond 
			if (dsq <= dist_to_bond*dist_to_bond) 
			{
				std::string pair = std::to_string(ip) + "-" + std::to_string(jp);
				if (bonded_pairs.count(pair) == 0)
				{
					if (particle_bond_type == "DEM") {
					    temp_bond = std::make_unique<Bond>(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ip, jp);
					} 
					else if (particle_bond_type == "linear") {
					    temp_bond = std::make_unique<Linear_elastic_bond>(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ip, jp);
					} 
					else {
					    throw std::invalid_argument("Wrong type of particle bond: " + particle_bond_type);
					}
					particle_bonds.push_back(std::move(temp_bond));

					bonded_pairs.insert(pair);
					// Increase the number of particle bonds
					particles.at(ip)->increase_bond_number();
					particles.at(jp)->increase_bond_number();			
				}
			}
 		}
	}
}

// Install new wall bonds or reinstall old one if within proximity limits
void CTC_sim::reinstall_wall_bonds(const double _k_n, const double _k_s, const double _sigma_max, 
            const double _tau_max, const double _R_bond, const double dist_to_bond)
{

 	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;
	// Coordinates of walls
	std::vector<double> walls = geom.get_box_limits();
	std::vector<double> y_coords{0,3};

	// Compute the distance between particles
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		for (size_t jp = 0; jp < 2; ++jp) 
		{
			// Compute the distance
			temp_particle_i = *(particles.at(ip));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);

			xj = temp_coords.at(0); 			
			yj = walls.at(y_coords.at(jp)); 
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);
	
			dsq = dx*dx + dy*dy;

			// If particle and wall are close enough and the bond is not present, form a bond 
			if (dsq <= dist_to_bond*dist_to_bond) 
			{
				std::string pair = std::to_string(ip) + "-wall";
				if (wall_bonded_pairs.count(pair) == 0)
				{
					// If this bond does not exist 
					wall_bonded_pairs.insert(pair);
					Wall_bond temp_bond(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ip, {xj, yj});
					wall_bonds.push_back(temp_bond);

					particles.at(ip)->set_wall_bond(true);
				}
			}
 		}
	}
}

void CTC_sim::initialize_bond_state(const double F_n0,
		const std::vector<double>& F_s0, const double M_b0)
{
	for (auto& bond_i : particle_bonds)
	{
		bond_i->set_bond_state(F_n0, F_s0, M_b0);
	}
}

void CTC_sim::initialize_wall_bond_state(const double F_n0,
		const std::vector<double>& F_s0, const double M_b0)
{
	for (auto& bond_i : wall_bonds)
	{
		bond_i.set_bond_state(F_n0, F_s0, M_b0);
	}
}

// Invoke boundary conditions
void CTC_sim::apply_boundary_conditions()
{	
	for (auto& particle : particles) {
		if (x_bc_type == "lj-walls"){ 
			// x condition (LJ)
			(boundaries.at(0))->apply_boundary_condition(particle->get_x_position(), particle->get_y_position());
			// y condition
			(boundaries.at(1))->apply_boundary_condition(particle->get_y_position());
		} else if (y_bc_type == "lj-walls") {
			// x condition
			(boundaries.at(0))->apply_boundary_condition(particle->get_x_position());
			// y condition
			(boundaries.at(1))->apply_boundary_condition(particle->get_y_position(), particle->get_x_position());
		} else {
			// x condition
			(boundaries.at(0))->apply_boundary_condition(particle->get_x_position());
			// y condition
			(boundaries.at(1))->apply_boundary_condition(particle->get_y_position());
		}
	} 
}

// Invoke xy boundary conditions (i.e. with two arguments)
void CTC_sim::apply_xy_boundary_conditions()
{	
	for (auto& particle : particles) {
		for (auto& boundary : xy_boundaries) {
			// x and y condition at the same time
			boundary->apply_boundary_condition(particle->get_x_position(), particle->get_y_position());
		}
	} 
}

// Invoke boundary conditions (mixed type)
void CTC_sim::apply_boundary_conditions_for_particle()
{	
	for (auto& particle : particles) {
		if (x_bc_type == "lj-walls"){ 
			// x condition (LJ)
			(boundaries.at(0))->apply_boundary_condition(particle->get_x_position(), particle->get_y_position());
			// y condition
			(boundaries.at(1))->apply_boundary_condition_for_particle(particle->get_y_position(), particle->get_wall_clearance());
		} else if (y_bc_type == "lj-walls") {
			// x condition
			(boundaries.at(0))->apply_boundary_condition_for_particle(particle->get_x_position(), particle->get_wall_clearance());
			// y condition
			(boundaries.at(1))->apply_boundary_condition(particle->get_y_position(), particle->get_x_position());
		} else {
			// x condition
			(boundaries.at(0))->apply_boundary_condition_for_particle(particle->get_x_position(), particle->get_wall_clearance());
			// y condition
			(boundaries.at(1))->apply_boundary_condition_for_particle(particle->get_y_position(), particle->get_wall_clearance());
		}
	} 
}
// Invoke xy boundary conditions (i.e. with two arguments) for mixed types
void CTC_sim::apply_xy_boundary_conditions_for_particle()
{	
	for (auto& particle : particles) {
		for (auto& boundary : xy_boundaries) {
			// x and y condition at the same time
			boundary->apply_boundary_condition_for_particle(particle->get_x_position(), particle->get_y_position(), particle->get_wall_clearance());
		}
	} 
}

// Apply external force field (displaces particles)
void CTC_sim::apply_external_force_field()
{
	double Fx = 0.0, Fy = 0.0, mass = 0.0;	
	for (auto& particle : particles) 
	{
		// Fetch particle mass
		mass = particle->get_mass();
		// Fetch force components for both coordinates
		// This needs to be improved, it assumes both forces depend only on y
		Fx = external_force.force_on_x(particle->get_y_position());
	   	Fy = external_force.force_on_y(particle->get_y_position());                         	
 		// Compute and apply displacements
		particle->correct_x_position(Fx*delta_t*delta_t/mass);
		particle->correct_y_position(Fy*delta_t*delta_t/mass);	
	}
}

// Find the average directions of the neighbors (particles within radius R)
void CTC_sim::compute_average_heading(const double R)
{
	// Variables for distance computation
	const double R_sq = R*R;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;
	
	// Computation of average heading
	double sin_theta = 0.0, cos_theta = 0.0;
	size_t neighbors = 0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i, temp_particle_j;

	// Find neighbors, compute average heading
	// Average includes oneself
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		neighbors = 0;
		sin_theta = 0.0;
		cos_theta = 0.0;
		for (size_t jp = 0; jp < particles.size(); ++jp) {	
			temp_particle_i = *(particles.at(ip));
			temp_particle_j = *(particles.at(jp));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);		
			dsq = dx*dx + dy*dy;

			// Neighbor found - count it and add to average heading
			if (dsq <= R_sq) {
				++neighbors;
				sin_theta += std::sin(temp_particle_j.get_heading());		
				cos_theta += std::cos(temp_particle_j.get_heading());
			}
 		}
		// Compute and save the new average heading
		particles.at(ip)->set_new_heading(std::atan2(sin_theta/neighbors, cos_theta/neighbors) + noise_mag*rng.get_random(-0.5, 0.5));
	}	

	// Update all headings with new values
	//std::for_each(particles.begin(), particles.end(), [](std::unique_ptr<Particle>& particle){ particle->update_heading(); });
}

// Find the average directions of the neighbors (particles within radius R) in the presence of an external flow
void CTC_sim::compute_average_heading(const double R, const bool needs_external_flow)
{
	// Variables for distance computation
	const double R_sq = R*R;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;
	
	// Computation of average heading
	double sin_theta = 0.0, cos_theta = 0.0, theta_j_hat = 0.0;
	double vel_neighbor = 0.0, theta_neighbor =0.0;
	double vx = 0.0, vy = 0.0;
	size_t neighbors = 0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i, temp_particle_j;

	// Find neighbors, compute average heading
	// Average includes oneself
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		neighbors = 0;
		sin_theta = 0.0;
		cos_theta = 0.0;
		for (size_t jp = 0; jp < particles.size(); ++jp) {	
			temp_particle_i = *(particles.at(ip));
			temp_particle_j = *(particles.at(jp));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);		
			dsq = dx*dx + dy*dy;

			// Neighbor found - count it and add to average heading
			if (dsq <= R_sq) {
				++neighbors;
				// Heading that accounts for external flow
				vel_neighbor = temp_particle_j.get_velocity_magnitude();
				theta_neighbor = temp_particle_j.get_heading();
				vx = external_flow.velocity_x(xj, yj);
				vy = external_flow.velocity_y(xj, yj);

				theta_j_hat = std::atan2(vel_neighbor*std::sin(theta_neighbor) + vy, 
											vel_neighbor*std::cos(theta_neighbor) + vx);

				sin_theta += std::sin(theta_j_hat);		
				cos_theta += std::cos(theta_j_hat);
			}
 		}
		// Compute and save the new average heading
		particles.at(ip)->set_new_heading(std::atan2(sin_theta/neighbors, cos_theta/neighbors) + noise_mag*rng.get_random(-0.5, 0.5));
	}	

	// Compute and store angular velocities
	// So lambda can capture it
	double dt = delta_t;
	std::for_each(particles.begin(), particles.end(), [&dt](std::unique_ptr<Particle>& particle){ particle->compute_angular_velocity(dt); });
}

void CTC_sim::compute_average_heading_two_types(const std::vector<double> R, const bool needs_external_flow)
{
	// Variables for distance computation
	const double RA = R.at(0);
	const double RB = R.at(1);
	const double R_sqA = RA*RA;
	const double R_sqB = RB*RB;
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;
	double R_sq = R_sqA;

	// Computation of average heading
	double sin_theta = 0.0, cos_theta = 0.0, theta_j_hat = 0.0;
	double vel_neighbor = 0.0, theta_neighbor =0.0;
	double vx = 0.0, vy = 0.0;
	size_t neighbors = 0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i, temp_particle_j;

	// Find neighbors, compute average heading
	// Average includes oneself
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		neighbors = 0;
		sin_theta = 0.0;
		cos_theta = 0.0;
		temp_particle_i = *(particles.at(ip));
		const size_t ti = temp_particle_i.get_type_id();
		if (ti == 0) {
			R_sq = R_sqA;
		} else {
			R_sq = R_sqB;
		}
		for (size_t jp = 0; jp < particles.size(); ++jp) {		
			temp_particle_j = *(particles.at(jp));
			
			// Same-type filter
			if (temp_particle_j.get_type_id() != ti) {
				continue;
			}
			// Positions
			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);		
			dsq = dx*dx + dy*dy;

			// Neighbor found - count it and add to average heading
			if (dsq <= R_sq) {
				++neighbors;
				// Heading that accounts for external flow
				vel_neighbor = temp_particle_j.get_velocity_magnitude();
				theta_neighbor = temp_particle_j.get_heading();
				
				if (needs_external_flow) {
					vx = external_flow.velocity_x(xj, yj);
					vy = external_flow.velocity_y(xj, yj);

					theta_j_hat = std::atan2(vel_neighbor*std::sin(theta_neighbor) + vy, 
												vel_neighbor*std::cos(theta_neighbor) + vx);
				} else {
					theta_j_hat = theta_neighbor;
				}
				
				sin_theta += std::sin(theta_j_hat);		
				cos_theta += std::cos(theta_j_hat);
			}
 		}
		if (ti == 0) {
			// Compute and save the new average heading for type A
			particles.at(ip)->set_new_heading(std::atan2(sin_theta/neighbors, cos_theta/neighbors) + noise_magA*rng.get_random(-0.5, 0.5));
		} else {
			// Compute and save the new average heading for type B
			particles.at(ip)->set_new_heading(std::atan2(sin_theta/neighbors, cos_theta/neighbors) + noise_magB*rng.get_random(-0.5, 0.5));
		}
	
	}	

	// Compute and store angular velocities
	// So lambda can capture it
	double dt = delta_t;
	std::for_each(particles.begin(), particles.end(), [&dt](std::unique_ptr<Particle>& particle){ particle->compute_angular_velocity(dt); });
	// Update all headings with new values
	std::for_each(particles.begin(), particles.end(), [](std::unique_ptr<Particle>& particle){ particle->update_heading(); });
}


// Move particles following their current velocity
void CTC_sim::propagate_particles()
{
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		particles.at(ip)->compute_velocity();
		particles.at(ip)->store_flow_contributions(delta_t);
	}
}

// Move particles following their current velocity and external flow
void CTC_sim::propagate_particles(const bool needs_external_flow)
{
	// Particle centroid
	double x = 0.0, y = 0.0;
	// Computation
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		x = particles.at(ip)->get_x_position();
		y = particles.at(ip)->get_y_position();

		particles.at(ip)->compute_velocity();

		particles.at(ip)->add_external_flow_contribution
			(external_flow.velocity_x(x, y), external_flow.velocity_y(x, y));

		particles.at(ip)->store_flow_contributions(delta_t);
	}
}

// Move away particles closer than d_min
void CTC_sim::remove_particle_overlaps(const double d_min) 
{ 
	// Variables for distance computation
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i, temp_particle_j;

	// Unit normal of the contact plane
	std::vector<double> nc{0.0, 0.0};

	// Compute the distance between particles
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		for (size_t jp = 0; jp < particles.size(); ++jp) {
		
			// Skip itself or if this pair was already considered	
			if (ip >= jp) {
				continue;
			}

			temp_particle_i = *(particles.at(ip));
			temp_particle_j = *(particles.at(jp));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);
	
			dsq = dx*dx + dy*dy;

			// Shift in opposite directions if closer than d_min 
			if (dsq < d_min*d_min) {
				basic_interactions.remove_overlaps(particles.at(ip), particles.at(jp), rng, dsq, d_min); 
			}
 		}
	}

	// Update all particle positions with new values
	std::for_each(particles.begin(), particles.end(), [](std::unique_ptr<Particle>& particle){ particle->update_position(); });

}

void CTC_sim::remove_particle_overlaps_two_types(const std::vector<double> d_min) 
{ 
	// Variables for distance computation
	const double d_AA = d_min.at(0);
	const double d_BB = d_min.at(1);
	double d_AB = 0.0;
	if (d_AA <= 0.0 || d_BB <= 0.0) {
		throw std::runtime_error("Error: Minimum distance for at least one particle type is not positive.");
	}
	// Average distance for different-type interactions
	if (d_min.size() > 2) {
		d_AB = d_min.at(2);
		if (d_AB <= 0.0) {
			throw std::runtime_error("Error: Minimum distance for different-type particles is not positive.");
		}
	} else {
		d_AB = 0.5*(d_AA + d_BB);
	}
	// Computation
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i, temp_particle_j;

	// Unit normal of the contact plane
	std::vector<double> nc{0.0, 0.0};

	// Compute the distance between particles
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		for (size_t jp = 0; jp < particles.size(); ++jp) {
		
			// Skip itself or if this pair was already considered	
			if (ip >= jp) {
				continue;
			}

			temp_particle_i = *(particles.at(ip));
			temp_particle_j = *(particles.at(jp));

			temp_coords = temp_particle_i.get_position();
			xi = temp_coords.at(0); 			
			yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			xj = temp_coords.at(0); 			
			yj = temp_coords.at(1);
	
			// Distance 
			dx = compute_x_distance(xi, xj);
			dy = compute_y_distance(yi, yj);
	
			dsq = dx*dx + dy*dy;

			// Shift in opposite directions if closer than threshold
			size_t ti = temp_particle_i.get_type_id();
			size_t tj = temp_particle_j.get_type_id();
			double d_target = d_AB;
			
			if (ti == tj) {
				if (ti == 0) {
					d_target = d_AA;
				} else {
					d_target = d_BB;
				}
			}

			if (dsq < d_target*d_target) {
				basic_interactions.remove_overlaps(particles.at(ip), particles.at(jp), rng, dsq, d_target); 
			}
 		}
	}

	// Update all particle positions with new values
	std::for_each(particles.begin(), particles.end(), [](std::unique_ptr<Particle>& particle){ particle->update_position(); });

}

// Compute bond forces and moments, check failure, correct particle properties
void CTC_sim::process_bonds()
{
	// Variables for distance computation
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i, temp_particle_j;
	std::vector<size_t> pIDs(2);

	std::vector<int> temp_failed;
	bool bond_failed = false;
	bool failed_in_shear = false, failed_in_tension = false;
	size_t failed_bonds = 0;

	// Iterate through all the bonds
	for (auto& bond_i : particle_bonds)
	{
		// Get the IDs of bonded particles
		pIDs = bond_i->get_particle_IDs();

		// Compute particle distance
		temp_particle_i = *(particles.at(pIDs.at(0)));
		temp_particle_j = *(particles.at(pIDs.at(1)));

		temp_coords = temp_particle_i.get_position();
		xi = temp_coords.at(0); 			
		yi = temp_coords.at(1);
		temp_coords = temp_particle_j.get_position();
		xj = temp_coords.at(0); 			
		yj = temp_coords.at(1);
	
		// Distance 
		dx = compute_x_distance(xi, xj);
		dy = compute_y_distance(yi, yj);
	
		dsq = dx*dx + dy*dy;

		// Break the bond if distance is too large
		if (std::sqrt(dsq) > max_bond_dist)
		{
			++failed_bonds;
			temp_failed.push_back(2);
			// Remove this pair
			std::string pair = std::to_string(pIDs.at(0)) + "-" + std::to_string(pIDs.at(1));
		 	size_t npairs = bonded_pairs.erase(pair);
			assert(npairs);
			// Reduce the number of particle bonds
			particles.at(pIDs.at(0))->decrease_bond_number();
			particles.at(pIDs.at(1))->decrease_bond_number();
			// Change bond status
			bond_i->set_failed(true);

			continue;			 
		}

		// Compute bond forces and moments
		bond_i->add_displacement_increments(particles.at(pIDs.at(0)), particles.at(pIDs.at(1)), rng, std::sqrt(dsq), delta_t);

		// Check for failure
		bond_failed = bond_i->bond_failed(failed_in_shear, failed_in_tension);
		if (bond_failed) {
			// If failed
			++failed_bonds;
			if (failed_in_shear) {
				temp_failed.push_back(1);
			} else {
				temp_failed.push_back(0);
			}
			// Remove this pair
			std::string pair = std::to_string(pIDs.at(0)) + "-" + std::to_string(pIDs.at(1));
		 	size_t npairs = bonded_pairs.erase(pair);
			assert(npairs);
			// Reduce the number of particle bonds
			particles.at(pIDs.at(0))->decrease_bond_number();
			particles.at(pIDs.at(1))->decrease_bond_number();			
		} else {
			// Still bonded - add bond contribution to particle properties
			bond_i->correct_particle_positions_and_headings(particles.at(pIDs.at(0)), particles.at(pIDs.at(1)), delta_t);		
		}
	}

	// Update if anything failed
	if (failed_bonds)
	{
		// Save the number of bonds that failed at this step
		failed_bond_numbers.push_back(failed_bonds);
		// Store the failure modes
		bond_failure_mode.push_back(temp_failed);	
		// Remove all with failed flag	
		particle_bonds.remove_if([](std::unique_ptr<Bond>& bond_i) {
    		return bond_i->failed_remove(); });
	}
}

// Compute wall bond forces and moments, check failure, correct particle properties
void CTC_sim::process_wall_bonds()
{
	// Variables for distance computation
	double xi = 0.0, xj = 0.0, yi = 0.0, yj =0.0; 
	double dx = 0.0, dy = 0.0, dsq = 0.0;

	// Temporary variables
	std::vector<double> temp_coords;
	Particle temp_particle_i;
	size_t pID = 0;

	std::vector<int> temp_failed;
	bool bond_failed = false;
	bool failed_in_shear = false, failed_in_tension = false;
	size_t failed_bonds = 0;

	// Iterate through all the bonds
	for (auto& bond_i : wall_bonds)
	{
		// Get the IDs of wall bonded particle
		pID = bond_i.get_particle_ID();

		// Compute particle distance
		temp_particle_i = *(particles.at(pID));

		temp_coords = temp_particle_i.get_position();
		xi = temp_coords.at(0); 			
		yi = temp_coords.at(1);

		xj = bond_i.get_x_coord(); 			
		yj = bond_i.get_y_coord();
	
		// Distance 
		dx = compute_x_distance(xi, xj);
		dy = compute_y_distance(yi, yj);
	
		dsq = dx*dx + dy*dy;

		// Break the bond if distance is too large
		if (std::sqrt(dsq) > max_wall_bond_dist)
		{
			++failed_bonds;
			temp_failed.push_back(2);

			// Remove this pair
			std::string pair = std::to_string(pID) + "-wall";
		 	size_t npairs = wall_bonded_pairs.erase(pair);
			assert(npairs);

			// Change bond status
			bond_i.set_failed(true);

			// Update particle status
			particles.at(pID)->set_wall_bond(false);

			continue;			 
		}

		// Compute bond forces and moments
		bond_i.add_displacement_increments(particles.at(pID), rng, std::sqrt(dsq), delta_t);

		// Check for failure
		bond_failed = bond_i.bond_failed(failed_in_shear, failed_in_tension);
		if (bond_failed) {
			// If failed
			++failed_bonds;
			if (failed_in_shear) {
				temp_failed.push_back(1);
			} else {
				temp_failed.push_back(0);
			}
			// Update particle status
			particles.at(pID)->set_wall_bond(false);
			// Remove this pair
			std::string pair = std::to_string(pID) + "-wall";
		 	size_t npairs = wall_bonded_pairs.erase(pair);
			assert(npairs);
		} else {
			// Still bonded - add bond contribution to particle properties
			bond_i.correct_particle_positions_and_headings(particles.at(pID), delta_t);		
		}
	}

	// Update if anything failed
	if (failed_bonds)
	{
		// Store the failure modes
		wall_bond_failure_mode.push_back(temp_failed);	
		// Remove all with failed flag	
		wall_bonds.remove_if([](Wall_bond bond_i){ return bond_i.failed_remove(); });
	}
}

// Compute and store lift force contribution for each particle
void CTC_sim::add_lift_force(const double Rp, const double rho, const double mu)
{
	double CL = 0.0, mass = 0.0, F_lift = 0.0;
	// Standard Saffman coefficient
    const double K_saffman = 6.46; 
    // Finite difference step size for calculating shear (1 micron is usually safe)
    const double dy_step = 1.0e-6;

	for (auto& particle : particles) 
	{
		// Fetch particle mass and radius
		mass = particle->get_mass();
		// Select CL based on current clustering status
		CL = static_cast<bool>(particle->get_bond_number()) ? CL_cluster : CL_single;

        // 1. Get Particle State
        double x = particle->get_x_position();
        double y = particle->get_y_position();
        double u_particle = particle->get_x_velocity();

        // 2. Get Fluid Velocity (Query the External Flow Object directly)
        double u_fluid = external_flow.velocity_x(x, y);

        // 3. Calculate Local Shear Rate (du/dy) numerically
        // We sample velocity slightly above and below to find the slope
        double u_upper = external_flow.velocity_x(x, y + dy_step);
        double u_lower = external_flow.velocity_x(x, y - dy_step);
        
        // Central Difference Approximation: slope = (f(x+h) - f(x-h)) / 2h
        double shear_rate = (u_upper - u_lower) / (2.0 * dy_step);
        double abs_shear = std::fabs(shear_rate);

        // 4. Calculate Slip Velocity
        double u_slip = u_fluid - u_particle;

        // 5. Compute Lift Force Magnitude
        // Formula: F ~ 6.46 * R^2 * sqrt(rho * mu * |shear|) * u_slip
        F_lift = CL * K_saffman * (Rp * Rp) * std::sqrt(rho * mu * abs_shear) * u_slip;

        // 6. Apply Direction
        // If shear is positive (velocity increases with Y), lagging particles go UP (+Y)
        // If shear is negative (velocity decreases with Y), lagging particles go DOWN (-Y)
        double direction_sign = (shear_rate > 0) ? 1.0 : -1.0;
        double force_y = F_lift * direction_sign;

		// Compute and store displacements
		particle->store_lift_contribution_x(0.0);
		particle->store_lift_contribution_y(force_y*delta_t*delta_t/mass);	
	}
}

// Mean velocity magnitude over all particles - not normalized
double CTC_sim::compute_average_velocity_magnitude()
{
	double vel_p = 0.0;
	const size_t Np = particles.size();
	for (size_t ip = 0; ip < Np; ++ip) {
		particles.at(ip)->compute_velocity();
		vel_p += (particles.at(ip)->get_velocity_magnitude());
	}	
	return vel_p/Np;
}

// Polarization 
double CTC_sim::compute_polarization()
{
    std::complex<double> total_phi(0.0, 0.0);
    const size_t Np = particles.size();
    for (size_t ip = 0; ip < Np; ++ip) {
        total_phi += std::polar(1.0, particles.at(ip)->get_heading());
    }
    return std::fabs(total_phi)/Np;
}

/// Compute the new positions and headings from all contributions
void CTC_sim::update_positions_and_headings()
{
	std::for_each(particles.begin(), particles.end(), [](std::unique_ptr<Particle>& particle)
	{ particle->update_heading(); particle->compute_new_position(); particle->reset_bond_contributions(); });
}

// -------- For now these two functions are here only for testing purposes

/*
// Install a bond between particle pair
void CTC_sim::create_a_bond(const size_t ID_1, const size_t ID_2, 
			const double _k_n, const double _k_s, const double _sigma_max,
            const double _tau_max, const double _R_bond)
{
	std::string pair = ID_1 < ID_2 ? std::to_string(ID_1) + "-" + std::to_string(ID_2) :
			std::to_string(ID_2) + "-" + std::to_string(ID_1); 
	bonded_pairs.insert(pair);
	
	Bond temp_bond(_k_n, _k_s, _sigma_max, _tau_max, _R_bond, ID_1, ID_2);
	particle_bonds.push_back(temp_bond); 
}

// Remove a bond between particle pair
// ----- This one doesn't actually remove the bond - just the pair information
void CTC_sim::remove_a_bond(const size_t ID_1, const size_t ID_2)
{
	std::string pair = 	ID_1 < ID_2 ? std::to_string(ID_1) + "-" + std::to_string(ID_2) :
			std::to_string(ID_2) + "-" + std::to_string(ID_1);
	size_t npairs = bonded_pairs.erase(pair);
	assert(npairs);
}
*/
