#include "../../../include/ctc_sim.h"

/*****************************************************
 * 
 * Simulation of a single CTC cluster in a 
 * microfluidic pipe (2D) with a 3x10^{-4} diameter
 *
 *****************************************************/

// Uses static variables because pointers to lambda don't capture
// https://stackoverflow.com/questions/28746744/passing-capturing-lambda-as-function-pointer
// See first (79+) comment (among others)
// To be refactored
static double flow_fluid_mu = 0.0, pipe_diameter = 0.0;
static double flow_pressure_drop = 0.0;

int main() {

	// Parameters
	LoadParameters ldp;
	
	// Load parameters that change during simulation	
	std::map <std::string, double> particle_parameters = ldp.load_parameter_map("parameters.txt");
	
	// Velocity field
	const bool needs_external_flow = true;

	// Overwrite the static variables
	pipe_diameter = particle_parameters.at("Pipe diameter");		
	flow_fluid_mu = particle_parameters.at("Viscosity of the fluid");	
	flow_pressure_drop = particle_parameters.at("Pressure drop");	

	// x component (uses static variables)
    auto vel_x = [](double x, double y) 
	{
   		return (flow_pressure_drop/(4*flow_fluid_mu)*
				(pipe_diameter*y-y*y));
    };
	
	// y component
	auto vel_y = [](double x, double y) { return 0.0; };

	// Check the flow		
	/*size_t N = 301;
	double dy = pipe_diameter / double (N-1);
	std::vector<double> u(N);
	for(size_t i=0; i<N; ++i){
		double y = i * dy;
		u.at(i) = vel_x(0.0, y); 
	}
	Easy_IO eflow_iovel("velocity_profile.txt");
    eflow_iovel.write_vector<double>(u);
	*/	

	// API
	CTC_sim eflow(vel_x, vel_y);

	// Time step
	double dt = particle_parameters.at("Time step");
	// Max simulation time
	size_t max_steps = static_cast<size_t>(particle_parameters.at("Maximum simulation time")/dt);

	// Noise and time step
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));

	// Geometry and boundary conditions
	std::string x_type{"periodic"}, y_type{"walls"};
	eflow.create_rectangular_enclosure(particle_parameters.at("Domain length x"), 
								particle_parameters.at("Pipe diameter"),
								particle_parameters.at("Domain center x"), 
								particle_parameters.at("Domain center y"));
	eflow.set_mixed_boundary_conditions(x_type, y_type, 
					particle_parameters.at("Minimum allowable distance with wall"));
	eflow.set_mixed_distances(x_type, "walls");

	// Particle initialization
	std::vector<double> all_lim{0.0, 0.0, particle_parameters.at("Initialization side"),
							    particle_parameters.at("Initialization side")};
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.set_initial_heading_interval(-EASY_PI/2.0, EASY_PI/2.0);
	
	// Initialize particles and bonds
	eflow.initialize_particles_from_limits(particle_parameters.at("Number of particles"), 
			x_lim, y_lim, particle_parameters.at("Velocity magnitude"), particle_parameters.at("Cell mass"));
	eflow.install_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));
	eflow.install_wall_bonds(particle_parameters.at("Wall bond normal stiffness"), 
						particle_parameters.at("Wall bond shear stiffness"),
						particle_parameters.at("Wall bond normal strength"),
						particle_parameters.at("Wall bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Wall bond max gap"));

	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);

	// Create bonds file
	std::string fbond = "bonds.txt";
	std::ofstream bonds_file(fbond);
	if (!bonds_file)
	{
		std::cerr << "Error: Could not open bonds.txt for writing." << std::endl;
		return 1;
	}

	// Create wall bonds file
	std::string fwbond = "wall_bonds.txt";
	std::ofstream wall_bonds_file(fwbond);
	if (!wall_bonds_file)
	{
		std::cerr << "Error: Could not open wall_bonds.txt for writing." << std::endl;
		return 1;
	}

	// Loop to run the repetitions
	for (size_t ir = 0; ir < static_cast<size_t>(particle_parameters.at("Number of repetitions")); ++ir) 
	{
		
		// Obtain valid set of particle pointers and bond references for this repetition
		const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
		// Bond pair data 
		const auto& bond_pairs = eflow.get_bonded_pairs_IDs();
		// Wall bond pair data 
		const auto& wall_bonded_pairs = eflow.get_wall_bonded_pairs_IDs();

		// Enforce boundary conditions
		eflow.apply_mixed_boundary_conditions();
		
		// Collect the initial x, y coordinates of particles
		for (size_t ip = 0; ip < particles.size(); ++ip) 
		{
			temp_info.at(0) = particles.at(ip)->get_x_position();
			temp_info.at(1) = particles.at(ip)->get_y_position();
			temp_info.at(2) = static_cast<double>(0);
			temp_info.at(3) = static_cast<double>(ir);
			particle_info.push_back(temp_info);
		}
		// Save bonded pairs IDs to file as lower ID-higher ID
		for (const auto& bond : bond_pairs) 
		{
			bonds_file << bond << " " << 0 << " " << ir << "\n";
		}
		// Save wall-bonded pairs IDs to file as cell ID-wall
		for (const auto& bond : wall_bonded_pairs) 
		{
			wall_bonds_file << bond << " " << 0 << " " << ir << "\n";
		}

		// Run the simulation from t = 1 to max_steps
		for (size_t it = 1; it <= max_steps; ++it) 
		{	
	
			// Align heading with neighboring cells
			eflow.compute_average_heading(particle_parameters.at("Interaction radius"),
					needs_external_flow);
	
			// Update particle positions
			eflow.propagate_particles(needs_external_flow);
	
			// Adjust positions to avoid overlap
			eflow.remove_particle_overlaps(particle_parameters.at("Minimum allowable particle distance"));

			// Compute bond forces, check for failure, add contributions
			eflow.process_bonds();			
	
			// Save active bond count
			eflow.track_bonds();

			// Compute bond forces, check for failure, add contributions
			eflow.process_wall_bonds();		

			// Enforce boundary conditions
			eflow.apply_mixed_boundary_conditions();

			// Install new bonds between cells and walls if any
			eflow.reinstall_wall_bonds(particle_parameters.at("Wall bond normal stiffness"), 
						particle_parameters.at("Wall bond shear stiffness"),
						particle_parameters.at("Wall bond normal strength"),
						particle_parameters.at("Wall bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Wall bond max gap"));
	
		}

		// Save final state
		for (size_t ip = 0; ip < particles.size(); ++ip) 
		{
			temp_info.at(0) = particles.at(ip)->get_x_position();
			temp_info.at(1) = particles.at(ip)->get_y_position();
			temp_info.at(2) = static_cast<double>(max_steps);
			temp_info.at(3) = static_cast<double>(ir);
			particle_info.push_back(temp_info);
		}
		// Save bonded pairs IDs to file as lower ID-higher ID
		for (const auto& bond : bond_pairs) 
		{
			bonds_file << bond << " " << max_steps << " " << ir << "\n";
		}
		// Save wall-bonded pairs IDs to file as cell ID-wall
		for (const auto& bond : wall_bonded_pairs) 
		{
			wall_bonds_file << bond << " " << max_steps << " " << ir << "\n";
		}

		// Re-initialize positions and heading angles
		eflow.reinitialize_particles(particle_parameters.at("Number of particles"), 
				x_lim, y_lim, particle_parameters.at("Velocity magnitude"), 
				particle_parameters.at("Cell mass"));
		
		// Reinstall bonds
		eflow.install_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));	
	
		eflow.install_wall_bonds(particle_parameters.at("Wall bond normal stiffness"), 
						particle_parameters.at("Wall bond shear stiffness"),
						particle_parameters.at("Wall bond normal strength"),
						particle_parameters.at("Wall bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Wall bond max gap"));
	}

	// Close bond files
	bonds_file.close();
	wall_bonds_file.close();

	// Save the results to file
	std::string fpos = "results.txt";
	Easy_IO eflow_io(fpos);
	eflow_io.write_vector<double>({particle_info});
	std::move(eflow_io);
	
}
