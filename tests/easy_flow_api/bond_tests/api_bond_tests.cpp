#include "../../../include/ctc_sim.h"
#include "../../../include/test_utils.h"
#include <set> 

/***************************************************** 
 *
 * Test suite for the ctc_sim API for bond usage 
 *	and formation. The results need to be inspected,
 * this test suite is not automated
 *
 *****************************************************/

// Uses static variables because pointers to lambda don't capture
// https://stackoverflow.com/questions/28746744/passing-capturing-lambda-as-function-pointer
// See first (79+) comment
// To be refactored
static double flow_fluid_mu = 0.0, flow_channel_width = 0.0;
static double flow_pressure_drop = 0.0;

int main() {

	// Parameters
	LoadParameters ldp;
	// Load parameters that change during simulation
	std::map <std::string, double> particle_parameters = ldp.load_parameter_map("parameters.txt");

	// Velocity field
	const bool needs_external_flow = true;
	// Overwrite the static variables
	flow_channel_width = particle_parameters.at("Domain length y");		
	flow_fluid_mu = particle_parameters.at("Viscosity of the fluid");
	flow_pressure_drop = particle_parameters.at("Pressure drop");
	// x component (uses static variables)
	auto vel_x = [](double x, double y) 
		{ return (((flow_channel_width*flow_channel_width/(2*flow_fluid_mu))
					*flow_pressure_drop)*(y/flow_channel_width)
					*(1-y/flow_channel_width)); }; 
	// y component
	auto vel_y = [](double x, double y) { return 0.0; };

	// Check the flow
	/*std::vector<double> y(100);
	std::iota(y.begin(), y.end(), 0.0);
	std::for_each(y.begin(), y.end(), [](double& yi) { yi *= flow_channel_width/100; });
	for (const auto& yi : y)
	{
		std::cout << yi << " " << vel_x(0.0, yi) << std::endl;
	}*/

	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Time stepping and collection settings
	// Time step
	double dt = particle_parameters.at("Time step");
	// Max simulation time
	size_t max_steps = static_cast<size_t>(particle_parameters.at("Maximum simulation time")/dt);

	// Collection times (time 0 is collected by default)
	std::vector<size_t> save_times(59);
	std::iota(save_times.begin(), save_times.end(), 1.0);
	// Convert to steps
	std::vector<size_t> save_steps_vector(save_times.size());
	std::transform(save_times.cbegin(), save_times.cend(), save_steps_vector.begin(), [&dt](size_t ti) 
		{ return static_cast<size_t>(static_cast<double>(ti)/dt); });
	std::set<size_t> save_steps(save_steps_vector.begin(), save_steps_vector.end());

	// Also noise
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));

	// Geometry and boundary conditions
	std::string x_type{"periodic"}, y_type{"walls"};
	eflow.create_rectangular_enclosure(particle_parameters.at("Domain length x"), 
								particle_parameters.at("Domain length y"),
								particle_parameters.at("Domain center x"), 
								particle_parameters.at("Domain center y"));
	eflow.set_mixed_boundary_conditions(x_type, y_type, 
					particle_parameters.at("Minimum allowable distance with wall"));
	eflow.set_mixed_distances(x_type, y_type);

	// Set parameters for particle initialization
	std::vector<double> all_lim{particle_parameters.at("Domain center x")
									- 0.6*2.0*particle_parameters.at("Particle radius"),
								particle_parameters.at("Domain center y")
									- 0.6*2.0*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center x")
									+ 0.6*2.0*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center y")
									+ 0.6*2.0*particle_parameters.at("Particle radius") };
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

	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);
	// For tracking test
	std::vector<size_t> current_bonds_with_time;

	// Outer loop - realizations 
	for (size_t ir = 0; ir < particle_parameters.at("Number of repetitions"); ++ir) {
		
		// Obtain valid set of particle pointers for this repetition
		const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();

		// Enforce boundary conditions
		eflow.apply_mixed_boundary_conditions();

		// Collect the initial x, y coordinates of particles 
		for (size_t ip = 0; ip < particles.size(); ++ip) {
			temp_info.at(0) = particles.at(ip)->get_x_position();
			temp_info.at(1) = particles.at(ip)->get_y_position();
			temp_info.at(2) = static_cast<double>(0);
			temp_info.at(3) = static_cast<double>(ir);
			particle_info.push_back(temp_info);
		}

		// Check if bonded pairs size is consistent with bond numbers
		const std::set<std::string>& bond_pairs = eflow.get_bonded_pairs_IDs();
		

		// Run the simulation from t = 1 to max_steps
		for (size_t it = 1; it < max_steps; ++it) {	
			
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

			// Compare to bond pairs
			current_bonds_with_time = eflow.get_bond_number_with_time();
			if (current_bonds_with_time.back() != bond_pairs.size())
			{
				std::cerr << "Number of bonds not matching number of ID pairs" << std::endl;
			}
			
			// Enforce boundary conditions
			eflow.apply_mixed_boundary_conditions();

			// Collect x, y coordinates of particles for specific time step in each realization
			if (save_steps.count(it) != 0) {
				for (size_t ip = 0; ip < particles.size(); ++ip) {
					temp_info.at(0) = particles.at(ip)->get_x_position();
					temp_info.at(1) = particles.at(ip)->get_y_position();
					temp_info.at(2) = static_cast<double>(it);
					temp_info.at(3) = static_cast<double>(ir);
					particle_info.push_back(temp_info);
				}
//				std::cout << "Saved time " << it*dt << std::endl; 
			}

		eflow.reinstall_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));
		}

		// Save active bond number to file
		std::vector<size_t> bonds_w_time = eflow.get_bond_number_with_time();
		Easy_IO eflow_io("active_bonds_w_time.txt");
		eflow_io.write_vector<size_t>({bonds_w_time});
	
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
	}

	// Save the results to file
	Easy_IO eflow_io("results.txt");
	eflow_io.write_vector<double>({particle_info});
	
}
