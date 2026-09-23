#include "../../include/ctc_sim.h"
#include "../../include/test_utils.h"
#include <set> 
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>

/*****************************************************
 * 
 * Simulation of CTC clusters in a microfluidic channel
 * with dimensions (0.12m x 0.001m x 0.003m) with 
 * 2D laminar flow in a rectangular section 
 * (0.12m /x/ x 0.01 /y/) as external flow and the 
 * bonded particle model as explicit 
 * bonding between cells.
 *
 * Source: 
 * doi: https://doi.org/10.1371/journal.pone.0245536
 *
 *****************************************************/

// Uses static variables because pointers to lambda don't capture
// https://stackoverflow.com/questions/28746744/passing-capturing-lambda-as-function-pointer
// See first (79+) comment (among others)
// To be refactored
static double flow_fluid_mu = 0.0, flow_channel_width = 0.0;
static double flow_pressure_drop = 0.0;
constexpr int N_MODES = 71;
static std::array<double, N_MODES> A, C;

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

	// Compute pressure drop
	double flow_channel_height = particle_parameters.at("Domain length z");
	double wss = particle_parameters.at("Wall shear stress");
	// Convert dyn to N/m2
	double wss_pa = wss * 0.1; 
	// Local correction factor for tau
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

	// Find maximum velocity		
	size_t N_grid = 301;
	double dy = flow_channel_width/static_cast<double>(N_grid-1), y = 0.0;
	std::vector<double> u_temp(N_grid);
	for(size_t i=0; i<N_grid; ++i){
		y = i * dy;
		u_temp[i] = vel_x(0.0, y); 
	}

	// External intialization (realization #1)
	int ini_result = 1;
	while (ini_result != 0) 
	{
	    ini_result = std::system("python3 make_initial.py");
	}
	// Load initial configuration 
	std::vector<double> x_ini, y_ini;
	std::ifstream file("initial_cells.txt");
	double xi, yi;
	while (file >> xi >> yi) {
	    x_ini.push_back(xi);
	    y_ini.push_back(yi);
	}
	file.close();
	
	// API
	CTC_sim eflow(vel_x, vel_y);

	// Time step
	double dt = particle_parameters.at("Time step");
	// Max simulation time
	size_t max_steps = static_cast<size_t>(particle_parameters.at("Maximum simulation time")/dt);
	// Data collection interval
	size_t collect_every = static_cast<size_t>(particle_parameters.at("Data collection interval")/dt);

	// Set up noise
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));

	// Set up the lift force coefficients
	eflow.set_lift_force_coefficient(particle_parameters.at("Single cell lift coefficient"),
				particle_parameters.at("Cell cluster lift coefficient"));

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
	const std::string bond_type = "linear";
	eflow.initialize_particles_from_list(particle_parameters.at("Number of particles"), 
			x_ini, y_ini, particle_parameters.at("Velocity magnitude"), 
			particle_parameters.at("Cell mass"));

	eflow.install_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"),
						particle_parameters.at("Max bond distance"), bond_type);
	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);

	// Create bonds file
	size_t wss_i = static_cast<size_t>(wss);
	std::string fbond = "bonds" + std::to_string(wss_i) + ".txt";
	std::ofstream bonds_file(fbond);
	if (!bonds_file){
		std::cerr << "Error: Could not open bonds.txt for writing." << std::endl;
		return 1;
	}

	// Accumulate bond failure codes across reps
	std::vector<int> all_bond_failures;

	// Outer loop - realizations 
	for (size_t ir = 0; ir < particle_parameters.at("Number of repetitions"); ++ir) {
		
		// Obtain valid set of particle pointers and bond references for this repetition
		const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
		// Bond pair data 
		const auto& bond_pairs = eflow.get_bonded_pairs_IDs();

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

			// Compute bond forces, check for failure, add contributions
			eflow.process_bonds();			

			// Lift force
			eflow.add_lift_force(particle_parameters.at("Particle radius"),
				particle_parameters.at("Density of the fluid"),
				particle_parameters.at("Viscosity of the fluid"));

			// Save active bond count
			eflow.track_bonds();

			// Update everything for this step
			eflow.update_positions_and_headings();

			// Adjust positions to avoid overlap
			eflow.remove_particle_overlaps(particle_parameters.at("Minimum allowable particle distance"));

			// Enforce boundary conditions
			eflow.apply_mixed_boundary_conditions();

			// Collect x, y coordinates of particles for specific time step in each realization
			if ((it % collect_every == 0) || (it == max_steps)) {
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

			// Reinstall bonds
			eflow.reinstall_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));
		}

		// Snapshot this realizations bond failure modes 
		for (const auto& batch : eflow.get_failure_mode_data())
		{
			for (int code : batch)
			{
				all_bond_failures.push_back(code);
			}
		}


		// External intialization (realizations > 1)
		ini_result = 1;
		while (ini_result != 0) 
		{
		    ini_result = std::system("python3 make_initial.py");
		}
		// Load initial configuration 
		x_ini.clear();
		y_ini.clear();
		file.clear();
		file.open("initial_cells.txt");
		while (file >> xi >> yi) 
		{
		    x_ini.push_back(xi);
		    y_ini.push_back(yi);
		}
		file.close();

		// Re-initialize positions and heading angles
		eflow.reinitialize_particles_from_list(particle_parameters.at("Number of particles"), 
				x_ini, y_ini, particle_parameters.at("Velocity magnitude"), 
				particle_parameters.at("Cell mass"));
		
		// Install bonds from scratch
		eflow.install_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"),
						particle_parameters.at("Max bond distance"), bond_type);
	}
	// Close bond file
	bonds_file.close();

	// Bond failure modes: one code per failure event (0=tension, 1=shear, 2=distance)
	std::string ffl = "bond_failures" + std::to_string(wss_i) + ".txt";
	std::ofstream wbf(ffl);
	for (int code : all_bond_failures){
		wbf << code << "\n";
	}
	wbf.close();

	// Save the results to file
	std::string fpos = "results" + std::to_string(wss_i) + ".txt";
	Easy_IO eflow_io(fpos);
	eflow_io.write_vector<double>({particle_info});
	
}	
