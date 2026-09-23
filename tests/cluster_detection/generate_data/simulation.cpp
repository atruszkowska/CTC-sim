#include "../../../include/ctc_sim.h"
#include "../../../include/test_utils.h"
#include <set> 
#include <cmath>

/***************************************************** 
 *
 * Verification of cluster detection algorithm  
 *
 *****************************************************/

// Uses static variables because pointers to lambda don't capture
// https://stackoverflow.com/questions/28746744/passing-capturing-lambda-as-function-pointer
// See first (79+) comment
static double flow_fluid_mu = 0.0, flow_channel_width = 0.0;
static double flow_pressure_drop = 0.0;
constexpr int N_MODES = 71;
static std::array<double, N_MODES> A, C;

int main() {

	// Parameters
	LoadParameters ldp;
	
	// Simulation parameters
	std::map <std::string, double> particle_parameters = ldp.load_parameter_map("parameters.txt");

   // Velocity field
   const bool needs_external_flow = true;

   // Overwrite the static variables
   flow_channel_width = particle_parameters.at("Domain length y");
   flow_fluid_mu = particle_parameters.at("Viscosity of the fluid");

   // Compute Pressure drop
   double flow_channel_height = particle_parameters.at("Domain length z");
   double wss = particle_parameters.at("Wall shear stress");
   double wss_pa = wss * 0.1; 
   double D = 0.3283; 
   flow_pressure_drop = (2 * wss_pa) / (D * flow_channel_height);

   // Precompute mode amplitudes & sech‑terms
   for(int m = 0; m < N_MODES; ++m) {
		// Odd mode index
		int n = 2*m + 1;
		double n3 = static_cast<double>(n*n*n);
		A.at(m) = (4.0 * flow_channel_width * flow_channel_width * flow_pressure_drop)
		        / (std::pow(EASY_PI, 3) * flow_fluid_mu * n3);
		C.at(m) = 1.0 / std::cosh(n * EASY_PI * flow_channel_height / (2.0 * flow_channel_width));
   }

   // x component (uses static variables)
   // mid‑plane (z=0) profile only
   auto vel_x = [](double x, double y) {
		double u = 0.0;
		for(int m = 0; m < N_MODES; ++m) {
		        int n = 2*m + 1;
		        u += A.at(m) * (1.0 - C.at(m)) * std::sin(n * EASY_PI * y / flow_channel_width);
		}
		return u;
   };
	
 	// y component
    auto vel_y = [](double x, double y) { return 0.0; };
	
	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Time stepping and collection settings
	// Time step
	double dt = particle_parameters.at("Time step");
	// Max simulation time
	size_t max_steps = static_cast<size_t>(particle_parameters.at("Maximum simulation time")/dt);

	// Collection times (time 0 is collected by default)
	std::vector<double> special_time_s = {
		0.6 * particle_parameters.at("Maximum simulation time"),
		particle_parameters.at("Maximum simulation time")
	};

	// Convert times in seconds to steps
	std::vector<size_t> save_steps_vector;
	save_steps_vector.reserve(special_time_s.size());

	for (double t_s : special_time_s) {
		// Convert time to number of steps
		size_t step_idx = static_cast<size_t>(std::floor(t_s / dt));
		save_steps_vector.push_back(step_idx);
	}

	std::set<size_t> save_steps(save_steps_vector.begin(), save_steps_vector.end());

	// Add noise
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
									- 0.5*particle_parameters.at("Initialization side"),
								0.0,
							    particle_parameters.at("Domain center x")
									+ 0.5*particle_parameters.at("Initialization side"),
							    particle_parameters.at("Domain length y")};
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.set_initial_heading_interval(-EASY_PI/2.0, EASY_PI/2.0);
	
	// Initialize particles and bonds
	eflow.initialize_particles_from_limits(particle_parameters.at("Number of particles"), 
			x_lim, y_lim, particle_parameters.at("Velocity magnitude"), 
			particle_parameters.at("Cell mass"));
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
	
	// Create bond file
	size_t wss_i = static_cast<size_t>(wss);
	std::string bnfile = "bonds" + std::to_string(wss_i) + ".txt";
	std::ofstream bonds_file(bnfile);
	if (!bonds_file) {
		std::cerr << "Error: Could not open bonds.txt for writing." << std::endl;
		return 1;
	} 	
	// Bond pair data 
	const auto& bond_pairs = eflow.get_bonded_pairs_IDs();

	// Sample initialization since each produces different number of clusters
	for (size_t ir = 0; ir < particle_parameters.at("Number of repetitions"); ++ir) {
	
		// Obtain a valid set of particle pointers
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
		// Save bonded pairs IDs to file as lower ID-higher ID
		for (const auto& bond : bond_pairs) {
			bonds_file << bond << " " << 0 << " " << ir << "\n";
		}
	
		// Run the simulation from t = 1 to max_steps
		for (size_t it = 1; it <= max_steps; ++it) {
			
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

             // Enforce boundary conditions
             eflow.apply_mixed_boundary_conditions();

			// Output bond number
			std::cout << bond_pairs.size() << std::endl;		

			// Collect x, y coordinates of particles for specific time step in each realization
			if (save_steps.count(it) != 0) {
				// Positions
				for (size_t ip = 0; ip < particles.size(); ++ip) {
					temp_info.at(0) = particles.at(ip)->get_x_position();
					temp_info.at(1) = particles.at(ip)->get_y_position();
					temp_info.at(2) = static_cast<double>(it);
					temp_info.at(3) = static_cast<double>(ir);
					particle_info.push_back(temp_info);
				}
				
				// Save bonded pairs IDs to file as lower ID-higher ID
				for (const auto& bond : bond_pairs) {
				    bonds_file << bond << " " << it << " " << ir << "\n";
				}
			}
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
        }


	// Save the results to file
	std::string posfile = "results" + std::to_string(wss_i) + ".txt";
	Easy_IO eflow_io(posfile);
	eflow_io.write_vector<double>({particle_info});
	
	// Close bond file
	bonds_file.close();	

}

