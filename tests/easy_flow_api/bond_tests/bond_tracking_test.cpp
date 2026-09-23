#include "../../../include/ctc_sim.h"
#include "../../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the ctc_sim API for tracking 
 *
 *****************************************************/

// Test 
bool bond_tracking_test();

int main()
{
 	test_pass(bond_tracking_test(), "Tracking bonds");
}

bool bond_tracking_test() 
{

	// Parameters
	LoadParameters ldp;
	// Load parameters that change during simulation
	std::map <std::string, double> particle_parameters = ldp.load_parameter_map("parameters.txt");

	// API
	CTC_sim eflow;
	
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

	// Create a bond between particles 1 and 3
	eflow.create_a_bond(1, 3, particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Bond max gap"));


	// Create a bond between particles 9 and 2
	eflow.create_a_bond(9, 2, particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Bond max gap"));


	// Verify correct bonds and IDs
	const std::set<std::string>& pairs_0 = eflow.get_bonded_pairs_IDs();
	if ((pairs_0.count("1-3") != 1) || (pairs_0.count("2-9") != 1))
	{
		for(auto it = pairs_0.cbegin(); it != pairs_0.cend(); it++)
	    {
	        std::cout << *it << std::endl;
	    }	

		std::cerr << "One of the bonds is missing" << std::endl;
		return false;
	}

	eflow.track_bonds();
	if (eflow.get_bond_number_with_time().at(0) != 2)
	{
		std::cerr << "There should be 2 bonds" << std::endl;
		return false;
	}
	
	// Remove the bond between particles 1 and 3
	eflow.remove_a_bond(3, 1);

	// Verify correct bonds and IDs
	if ((pairs_0.count("1-3") != 0))
	{
		for(auto it = pairs_0.cbegin(); it != pairs_0.cend(); it++)
	    {
	        std::cout << *it << std::endl;
	    }	

		std::cerr << "The 1-3 bond is not removed" << std::endl;
		return false;
	}

	// Remove the bond between particles 9 and 2
	eflow.remove_a_bond(2, 9);

	// Verify correct bonds and IDs
	if ((pairs_0.count("2-9") != 0))
	{
		for(auto it = pairs_0.cbegin(); it != pairs_0.cend(); it++)
	    {
	        std::cout << *it << std::endl;
	    }	

		std::cerr << "The 2-9 bond is not removed" << std::endl;
		return false;
	}

	return true;
}
