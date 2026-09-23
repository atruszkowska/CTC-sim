#include "../../include/ctc_sim.h"
#include "../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the ctc_sim API
 *
 *****************************************************/

static double flow_test_mu = 0.0, flow_test_dpdL = 0.0;
static double flow_test_channel_radius_squared =0.0, flow_test_channel_radius = 0.0;
static double F_x_max =0.0, F_y_max = 0.0, step_size = 0.0,	domain_height = 0.0;
static double flow_fluid_mu = 0.0, flow_channel_width = 0.0;
static double flow_pressure_drop = 0.0;

// Tests

bool geometry_tests();
bool initialization_tests();
bool re_initialization_tests();
bool initialization_tests_lists();
bool initialization_tests_round();
bool re_initialization_tests_round();
bool initialization_tests_two_types();
bool re_initialization_tests_two_types();
bool initialization_tests_round_two_types();
bool re_initialization_tests_round_two_types();
bool boundary_conditions_tests();
bool neighbor_heading_tests();
bool particle_propagation_tests();
bool particle_interactions_tests();
//bool test_wrong_boundary_types();
bool test_wrong_distance_types();
bool external_flow();
bool external_force();
bool lift_force_test();
bool lift_force_contribution_test();
bool all_contributions_test();
//bool dist_bond_breaking_test();
bool dist_wall_bond_breaking_test();

// Supporting functions
// Makes the rectangular domain and compares the resulting object with expectation
bool make_and_check_rectangle(const std::vector<double>&, const std::vector<double>&, 
					const std::vector<double>&, const std::vector<double>&);
// Makes and tests round enclosure
bool make_and_check_circle(const double, const std::vector<double>&, const double, const std::vector<double>&);
// Setup and test different booundary conditions
bool make_and_test_boundary_conditions(const std::string&, const std::string&);
bool make_and_test_boundary_conditions_lj(const std::string&, const std::string&);
bool make_and_test_xy_boundary_conditions();
// Sets up and tests heading computations for different boundary types
bool make_and_test_heading(const std::string&, const std::string&);
bool make_and_test_heading_two_types(const std::string&, const std::string&);
// Sets up and tests interactions computations for different boundary types
bool make_and_test_interactions(const std::string&, const std::string&);
bool make_and_test_interactions_two_types(const std::string&, const std::string&);
// Returns the number of overlapping particle pairs 
int count_overlaps(CTC_sim& eflow, const double d_min);
int count_overlaps(CTC_sim& eflow, const std::vector<double> d_mins);
// Returns distance between particle pairs
double compute_particle_distance(const Particle& pt1, const Particle& pt2);
// Simulate particles with external flow in a 2D channel
// Arguments are particle and external flow velocity magnitudes (v and u) and noise (eta)
// Returns particle coordinates and velocities
std::vector<std::vector<double>> apply_external_flow(const double v, const double u, const double eta);

int main()
{
	test_pass(geometry_tests(), "Geometry API");
 	test_pass(initialization_tests(), "Initialize particles");
	test_pass(re_initialization_tests(), "Reinitialize particles");
 	test_pass(initialization_tests_lists(), "Initialize and reinitialize particles from a list");
 	test_pass(initialization_tests_round(), "Initialize particles in a circular enclosure");
	test_pass(re_initialization_tests_round(), "Reinitialize particles in a circular enclosure");
	test_pass(initialization_tests_two_types(), "Initialize particles of two different types");
	test_pass(re_initialization_tests_two_types(), "Re-initialize particles of two different types");
	test_pass(initialization_tests_round_two_types(), "Initialize particle types in a circular enclosure");
	test_pass(re_initialization_tests_round_two_types(), "Re-initialize particle types in a circular enclosure");
	test_pass(boundary_conditions_tests(), "Set up and use boundary conditions");
	test_pass(neighbor_heading_tests(), "Computation of average heading based on the neighbors");
	test_pass(particle_propagation_tests(), "Particle propagation");
	test_pass(particle_interactions_tests(), "Particle interactions - removal of overlaps");
//	test_pass(test_wrong_boundary_types(), "Wrong types of boundary conditions");
	test_pass(test_wrong_distance_types(), "Wrong types of distance computation approaches");
	test_pass(external_flow(), "Particles in an external flow");
	test_pass(external_force(), "Particles in an external force field");
	test_pass(lift_force_test(), "Lift force and new computation reference");
	test_pass(lift_force_contribution_test(), "Lift force contribution");
	test_pass(all_contributions_test(), "Lift force, bonds, and self propelling contributions");
//	test_pass(dist_bond_breaking_test(), "Long bonds breaking");
//	test_pass(dist_wall_bond_breaking_test(), "Long wall bonds breaking");
}

// Construct a box and compare with expected
// There are three box and three well types considered
bool geometry_tests()
{
	// Input data
	const std::vector<std::vector<double>> all_dimensions{{10.2, 10.2}, 
										{12.8, 21.0}, {8.1, 2.42}};		
	const std::vector<std::vector<double>> all_centers{{0.1, 0.1}, 
										{6.4, 10.5}, {-1.0, 0.0}};

	const std::vector<std::vector<double>> expected_dimensions{{10.2, 10.2}, 
										{12.8, 21.0}, {8.1, 2.42}};

	const std::vector<std::vector<double>> expected_bounds{{-5, -5, 5.2, 5.2}, 
										{0, 0, 12.8, 21}, {-5.05, -1.21, 3.05, 1.21}};

	// Run a test on each set
	for (size_t i = 0; i < all_dimensions.size(); ++i) {
		if (!make_and_check_rectangle(all_dimensions.at(i), all_centers.at(i),
			expected_dimensions.at(i), expected_bounds.at(i))) {
			std::cerr << "Box making test failed for dataset " << i << std::endl;
			return false;
		}
	}

	// Input data - circle
	const std::vector<double> all_dimensions_circle{10.2, 2.1e-7, 1.6};		
	const std::vector<std::vector<double>> all_centers_circle{{0.1, 0.1}, 
										{6.4, 10.5}, {-1.0, 0.0}};
	const std::vector<double> expected_dimensions_circle{10.2, 2.1e-7, 1.6};		
	const std::vector<std::vector<double>> expected_centers_circle{{0.1, 0.1}, 
										{6.4, 10.5}, {-1.0, 0.0}};

	// Run a test on each set
	for (size_t i = 0; i < all_dimensions_circle.size(); ++i) {
		if (!make_and_check_circle(all_dimensions_circle.at(i), all_centers_circle.at(i),
			expected_dimensions_circle.at(i), expected_centers_circle.at(i))) {
			std::cerr << "Circle making test failed for dataset " << i << std::endl;
			return false;
		}
	}

	return true;
}

bool make_and_check_rectangle(const std::vector<double>& dimensions,
					const std::vector<double>& center,	
					const std::vector<double>& expected_dimensions, 
					const std::vector<double>& expected_bounds)
{
	double tol = 1e-5;
	CTC_sim eflow;

	eflow.create_rectangular_enclosure(dimensions.at(0), dimensions.at(1), center.at(0), center.at(1));

	// Bounds check
	const std::vector<double> bounds = eflow.get_enclosure_limits();
	if (!std::equal(bounds.cbegin(), bounds.cend(), expected_bounds.cbegin(), 
			[tol](double comp, double exp){ return float_equality(comp, exp, tol); })) {
		std::cerr << "Bounds not equal to expected"<< std::endl;
		return false;
	}
	
	return true;
}

bool make_and_check_circle(const double radius,
					const std::vector<double>& center,	
					const double expected_radius, 
					const std::vector<double>& expected_center)
{
	double tol = 1e-5;
	CTC_sim eflow;

	eflow.create_circular_enclosure(radius, center.at(0), center.at(1));

	// Dimensions and center check
	if (!(float_equality(eflow.get_enclosure_radius(), expected_radius, tol))) {
		std::cerr << "Circle radius does not match expected" << std::endl;
		return false;
	}
	const std::vector<double> cur_center = eflow.get_enclosure_center();
	if (!std::equal(cur_center.cbegin(), cur_center.cend(), expected_center.cbegin(), 
			[tol](double comp, double exp){ return float_equality(comp, exp, tol); })) {
		std::cerr << "Center coordinates not equal to expected"<< std::endl;
		return false;
	}
	
	return true;
}

// Tests initialization of particles with the same speed
bool initialization_tests()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	const double tol = 0.1; 	

	// Create the geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	// Average values following uniform distribution
	const double av_x_pos = 0.5*(x_lim.at(0) + x_lim.at(1));
	const double av_y_pos = 0.5*(y_lim.at(0) + y_lim.at(1));
	const double av_heading = 0.5*(-EASY_PI/2.0 + EASY_PI/2.0);
	
	// Initialization settings
	const size_t Np = 100000;
	const double v_mag = 0.1;

	// Reusable heading
	double theta = 0.0;

	// Initialize 1000 particles
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag);

	// Collect positions and headings
	std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		all_x_coords.push_back(coords.at(0));
		all_y_coords.push_back(coords.at(1));
		theta = (particles.at(ip))->get_heading();
		all_headings.push_back(theta);
		// Check limits
		if ((coords.at(0) < x_lim.at(0)) ||
			(coords.at(0) > x_lim.at(1)) ||
			(coords.at(1) < y_lim.at(0)) ||
			(coords.at(1) > y_lim.at(1))) {
				std::cerr << "x or y coordinates of particle " 
						  << ip << " are out of bounds" << std::endl;
				return false;
		}
		// Check particle velocity
 		vel = (particles.at(ip))->get_velocity();
		if (!(float_equality(vel.at(0), v_mag*std::cos(theta), tol))) {
			std::cerr << "x velocity component does not match expected" << std::endl;
			return false;
		}
		if (!(float_equality(vel.at(1), v_mag*std::sin(theta), tol))) {
			std::cerr << "y velocity component does not match expected" << std::endl;
			return false;
		}
	}

	// Check if the values follow uniform distribution
	const double ptcl_x_pos = std::accumulate(all_x_coords.begin(), all_x_coords.end(), 0.0)/Np;
	const double ptcl_y_pos = std::accumulate(all_y_coords.begin(), all_y_coords.end(), 0.0)/Np;
	const double ptcl_heading = std::accumulate(all_headings.begin(), all_headings.end(), 0.0)/Np;
	if (!(float_equality(ptcl_x_pos, av_x_pos, tol))) {
		std::cerr << "Average x coordinate " << ptcl_x_pos  
				  << " does not match expected " <<  av_x_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_y_pos, av_y_pos, tol))) {
		std::cerr << "Average y coordinate " << ptcl_y_pos  
				  << " does not match expected " <<  av_y_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_heading, av_heading, tol))) {
		std::cerr << "Average heading " << ptcl_heading  
				  << " does not match expected " <<  av_heading << std::endl;
		return false;
	}

	return true;	
	 
}

// Tests re-initialization of particles with the same speed
bool re_initialization_tests()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	const double tol = 0.1; 	

	// Create the geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	// Average values following uniform distribution
	const double av_x_pos = 0.5*(x_lim.at(0) + x_lim.at(1));
	const double av_y_pos = 0.5*(y_lim.at(0) + y_lim.at(1));
	const double av_heading = 0.5*(-EASY_PI/2.0 + EASY_PI/2.0);
	
	// Initialization settings
	const size_t Np = 100000;
	const double v_mag = 0.1;

	// Reusable heading
	double theta = 0.0;

	// Initialize 1000 particles
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag);
	eflow.reinitialize_particles(Np, x_lim, y_lim, v_mag);

	// Collect positions and headings
	std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		all_x_coords.push_back(coords.at(0));
		all_y_coords.push_back(coords.at(1));
		theta = (particles.at(ip))->get_heading();
		all_headings.push_back(theta);
		// Check limits
		if ((coords.at(0) < x_lim.at(0)) ||
			(coords.at(0) > x_lim.at(1)) ||
			(coords.at(1) < y_lim.at(0)) ||
			(coords.at(1) > y_lim.at(1))) {
				std::cerr << "x or y coordinates of particle " 
						  << ip << " are out of bounds" << std::endl;
				return false;
		}
		// Check particle velocity
 		vel = (particles.at(ip))->get_velocity();
		if (!(float_equality(vel.at(0), v_mag*std::cos(theta), tol))) {
			std::cerr << "x velocity component does not match expected" << std::endl;
			return false;
		}
		if (!(float_equality(vel.at(1), v_mag*std::sin(theta), tol))) {
			std::cerr << "y velocity component does not match expected" << std::endl;
			return false;
		}
	}

	// Check if the values follow uniform distribution
	const double ptcl_x_pos = std::accumulate(all_x_coords.begin(), all_x_coords.end(), 0.0)/Np;
	const double ptcl_y_pos = std::accumulate(all_y_coords.begin(), all_y_coords.end(), 0.0)/Np;
	const double ptcl_heading = std::accumulate(all_headings.begin(), all_headings.end(), 0.0)/Np;
	if (!(float_equality(ptcl_x_pos, av_x_pos, tol))) {
		std::cerr << "Average x coordinate " << ptcl_x_pos  
				  << " does not match expected " <<  av_x_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_y_pos, av_y_pos, tol))) {
		std::cerr << "Average y coordinate " << ptcl_y_pos  
				  << " does not match expected " <<  av_y_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_heading, av_heading, tol))) {
		std::cerr << "Average heading " << ptcl_heading  
				  << " does not match expected " <<  av_heading << std::endl;
		return false;
	}

	return true; 
}

bool initialization_tests_lists()
{
    CTC_sim eflow;

    const size_t Np = 3;
    const std::vector<double> x_pos = {1.0, 2.0, 3.0};
    const std::vector<double> y_pos = {4.0, 5.0, 6.0};
    const double vmag = 2.5;
    const double mass = 1.7;

    // Test initialization
    eflow.initialize_particles_from_list(
        Np, x_pos, y_pos, vmag, mass
    );

    const auto& particles = eflow.get_particles();

    if (particles.size() != Np)
        return false;

    for (size_t i = 0; i < Np; ++i)
    {
        const std::vector<double> position = particles.at(i)->get_position();

        if (position.at(0) != x_pos.at(i))
            return false;

        if (position.at(1) != y_pos.at(i))
            return false;

        if (particles.at(i)->get_velocity_magnitude() != vmag)
            return false;

        if (particles.at(i)->get_mass() != mass)
            return false;
    }

    // Test reinitialization
    const size_t Np_new = 2;
    const std::vector<double> x_new = {10.0, 20.0};
    const std::vector<double> y_new = {30.0, 40.0};
    const double vmag_new = 5.0;
    const double mass_new = 3.0;

    eflow.reinitialize_particles_from_list(
        Np_new, x_new, y_new, vmag_new, mass_new
    );

    const auto& particles_new = eflow.get_particles();

    // Verify that the old particles were replaced
    if (particles_new.size() != Np_new)
        return false;

    for (size_t i = 0; i < Np_new; ++i)
    {
        const std::vector<double> position = particles_new.at(i)->get_position();

        if (position.at(0) != x_new.at(i))
            return false;

        if (position.at(1) != y_new.at(i))
            return false;

        if (particles_new.at(i)->get_velocity_magnitude() != vmag_new)
            return false;

        if (particles_new.at(i)->get_mass() != mass_new)
            return false;
    }

    return true;
}

// Tests initialization of particles with the same speed in a circular enclosure
bool initialization_tests_round()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	double tol = 1e-5; 	

	// Geometry
	const double R_enclosure = 2.1;
	std::vector<double> enclosure_center = {0.5, 0.1};
	const double av_x_pos = 0.5, av_y_pos = 0.1, av_heading = 0.0; 

	// Initialization settings
	const size_t Np = 100000;
	const double v_mag = 0.1;

	// Reusable particle information
	double x0 = 0.0, y0 = 0.0, theta = 0.0;

	// Initialize particles
	eflow.initialize_particles_round_enclosure(Np, enclosure_center, R_enclosure, v_mag);

	// Check and collect positions and headings
	std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		x0 = coords.at(0);
		y0 = coords.at(1);
		all_x_coords.push_back(x0);
		all_y_coords.push_back(y0);
		theta = (particles.at(ip))->get_heading();
		all_headings.push_back(theta);
		// Check limits
		x0 -= enclosure_center.at(0);
		y0 -= enclosure_center.at(1);
		if ((x0*x0 + y0*y0) > R_enclosure*R_enclosure) {
				std::cerr << "x or y coordinates of particle " 
						  << ip << " are out of bounds" << std::endl;
				return false;
		}
		// Check particle velocity
 		vel = (particles.at(ip))->get_velocity();
		if (!(float_equality(vel.at(0), v_mag*std::cos(theta), tol))) {
			std::cerr << "x velocity component does not match expected" << std::endl;
			return false;
		}
		if (!(float_equality(vel.at(1), v_mag*std::sin(theta), tol))) {
			std::cerr << "y velocity component does not match expected" << std::endl;
			return false;
		}
	}

	// Check if the values follow uniform distribution
	const double ptcl_x_pos = std::accumulate(all_x_coords.begin(), all_x_coords.end(), 0.0)/Np;
	const double ptcl_y_pos = std::accumulate(all_y_coords.begin(), all_y_coords.end(), 0.0)/Np;
	const double ptcl_heading = std::accumulate(all_headings.begin(), all_headings.end(), 0.0)/Np;
	tol = 0.1;
	if (!(float_equality(ptcl_x_pos, av_x_pos, tol))) {
		std::cerr << "Average x coordinate " << ptcl_x_pos  
				  << " does not match expected " <<  av_x_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_y_pos, av_y_pos, tol))) {
		std::cerr << "Average y coordinate " << ptcl_y_pos  
				  << " does not match expected " <<  av_y_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_heading, av_heading, tol))) {
		std::cerr << "Average heading " << ptcl_heading  
				  << " does not match expected " <<  av_heading << std::endl;
		return false;
	}

	return true;	
	 
}

// Tests re-initialization of particles with the same speed in a round enclosure
bool re_initialization_tests_round()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	double tol = 1e-5; 	

	// Geometry
	const double R_enclosure = 2.1;
	std::vector<double> enclosure_center = {0.5, 0.1};
	const double av_x_pos = 0.5, av_y_pos = 0.1, av_heading = 0.0; 

	// Initialization settings
	const size_t Np = 100000;
	const double v_mag = 0.1;

	// Reusable particle information
	double x0 = 0.0, y0 = 0.0, theta = 0.0;

	// Initialize particles
	eflow.initialize_particles_round_enclosure(Np, enclosure_center, R_enclosure, v_mag);
	eflow.reinitialize_particles(Np, enclosure_center, R_enclosure, v_mag);

	// Check and collect positions and headings
	std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		x0 = coords.at(0);
		y0 = coords.at(1);
		all_x_coords.push_back(x0);
		all_y_coords.push_back(y0);
		theta = (particles.at(ip))->get_heading();
		all_headings.push_back(theta);
		// Check limits
		x0 -= enclosure_center.at(0);
		y0 -= enclosure_center.at(1);
		if ((x0*x0 + y0*y0) > R_enclosure*R_enclosure) {
				std::cerr << "x or y coordinates of particle " 
						  << ip << " are out of bounds" << std::endl;
				return false;
		}
		// Check particle velocity
 		vel = (particles.at(ip))->get_velocity();
		if (!(float_equality(vel.at(0), v_mag*std::cos(theta), tol))) {
			std::cerr << "x velocity component does not match expected" << std::endl;
			return false;
		}
		if (!(float_equality(vel.at(1), v_mag*std::sin(theta), tol))) {
			std::cerr << "y velocity component does not match expected" << std::endl;
			return false;
		}
	}

	// Check if the values follow uniform distribution
	const double ptcl_x_pos = std::accumulate(all_x_coords.begin(), all_x_coords.end(), 0.0)/Np;
	const double ptcl_y_pos = std::accumulate(all_y_coords.begin(), all_y_coords.end(), 0.0)/Np;
	const double ptcl_heading = std::accumulate(all_headings.begin(), all_headings.end(), 0.0)/Np;
	tol = 0.1;
	if (!(float_equality(ptcl_x_pos, av_x_pos, tol))) {
		std::cerr << "Average x coordinate " << ptcl_x_pos  
				  << " does not match expected " <<  av_x_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_y_pos, av_y_pos, tol))) {
		std::cerr << "Average y coordinate " << ptcl_y_pos  
				  << " does not match expected " <<  av_y_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_heading, av_heading, tol))) {
		std::cerr << "Average heading " << ptcl_heading  
				  << " does not match expected " <<  av_heading << std::endl;
		return false;
	}

	return true;
}

bool initialization_tests_two_types()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	const double tol = 0.1; 	

	// Create the geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	// Average values following uniform distribution
	//const double av_x_pos = 0.5*(x_lim.at(0) + x_lim.at(1));
	//const double av_y_pos = 0.5*(y_lim.at(0) + y_lim.at(1));
	//const double av_heading = 0.5*(-EASY_PI/2.0 + EASY_PI/2.0);
	
	// Initialization settings
	const std::vector<size_t> Np_types = {1000, 800};
	const std::vector<double> v_mags = {0.1, 1.0};
	const std::vector<double> R_mins = {0.05, 0.2};
	const std::vector<size_t> type_ids = {0, 1};
	const size_t Np = std::accumulate(Np_types.begin(), Np_types.end(), 0);
	size_t countA = 0;
	size_t countB = 0;

	// Reusable heading
	double theta = 0.0;

	// Initialize particles
	eflow.initialize_particles_from_limits(Np_types, x_lim, y_lim, v_mags, type_ids, R_mins);
	
	// Collect particles
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();	
	// Check particle type count	
	for (size_t ip = 0; ip < particles.size(); ++ ip) {
		if (particles.at(ip)->is_type(0)){
			countA += 1;
		} else if (particles.at(ip)->is_type(1)){
			countB += 1;
		} else {
			std::cerr << "Error wrong particle type found" << std::endl;
			return false;
		}
	}
	
	std::cout << "Count A: " << countA << " and expected: " << Np_types.at(0) << std::endl;
	std::cout << "Count B: " << countB << " and expected: " << Np_types.at(1) << std::endl;
	if (!(countA == Np_types.at(0) && countB == Np_types.at(1))) {
		std::cerr << "Error with assigning particle type" << std::endl;
		return false;
	}


	// Collect positions and headings
	std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		all_x_coords.push_back(coords.at(0));
		all_y_coords.push_back(coords.at(1));
		theta = (particles.at(ip))->get_heading();
		all_headings.push_back(theta);
		// Check limits
		if ((coords.at(0) < x_lim.at(0)) ||
			(coords.at(0) > x_lim.at(1)) ||
			(coords.at(1) < y_lim.at(0)) ||
			(coords.at(1) > y_lim.at(1))) {
				std::cerr << "x or y coordinates of particle " 
						  << ip << " are out of bounds" << std::endl;
				return false;
		}
		// Check particle velocity
		auto type = (particles.at(ip))->get_type_id();
 		vel = (particles.at(ip))->get_velocity();
		if (type == type_ids.at(0)) {
			if (!(float_equality(vel.at(0), v_mags.at(0)*std::cos(theta), tol))) {
				std::cerr << "x velocity component of type A does not match expected" << std::endl;
				return false;
			}
			if (!(float_equality(vel.at(1), v_mags.at(0)*std::sin(theta), tol))) {
				std::cerr << "y velocity component of type A does not match expected" << std::endl;
				return false;
			}
		} else if (type == type_ids.at(1)) {
			if (!(float_equality(vel.at(0), v_mags.at(1)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type B does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(1)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type B does not match expected" << std::endl;
                                return false;
                        }	
		
		} else {
			std::cerr << "Wrong particle type" << std::endl;
		}
	}
	
	return true;	
}

// Tests initialization of particles with the same speed in a circular enclosure
bool initialization_tests_round_two_types()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	double tol = 1e-5; 	

	// Geometry
	const double R_enclosure = 2.1;
	std::vector<double> enclosure_center = {0.5, 0.1};
	const double av_x_pos = 0.5, av_y_pos = 0.1, av_heading = 0.0; 

	// Initialization settings
	const std::vector<size_t> Np_types = {1000, 800};
        const std::vector<double> v_mags = {0.1, 1.0};
        const std::vector<double> R_mins = {0.05, 0.2};
        const std::vector<size_t> type_ids = {0, 1};
        const size_t Np = std::accumulate(Np_types.begin(), Np_types.end(), 0);
        size_t countA = 0;
        size_t countB = 0;

        // Reusable particle information
	double x0 = 0.0, y0 = 0.0, theta = 0.0;

	// Initialize particles
	eflow.initialize_particles_round_enclosure(Np_types, enclosure_center, R_enclosure, v_mags, type_ids, R_mins);

	// Check particle type count    
        const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < particles.size(); ++ ip) {
                if (particles.at(ip)->is_type(0)){
                        countA += 1;
                } else if (particles.at(ip)->is_type(1)){
                        countB += 1;
                } else {
                        std::cerr << "Error wrong particle type found" << std::endl;
                        return false;
                }
        }

        std::cout << "Count A: " << countA << " and expected: " << Np_types.at(0) << std::endl;
        std::cout << "Count B: " << countB << " and expected: " << Np_types.at(1) << std::endl;
        if (!(countA == Np_types.at(0) && countB == Np_types.at(1))) {
                std::cerr << "Error with assigning particle type" << std::endl;
                return false;
        }


	// Check and collect positions and headings
	std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		x0 = coords.at(0);
		y0 = coords.at(1);
		all_x_coords.push_back(x0);
		all_y_coords.push_back(y0);
		theta = (particles.at(ip))->get_heading();
		all_headings.push_back(theta);
		// Check limits
		x0 -= enclosure_center.at(0);
		y0 -= enclosure_center.at(1);
		if ((x0*x0 + y0*y0) > R_enclosure*R_enclosure) {
				std::cerr << "x or y coordinates of particle " 
						  << ip << " are out of bounds" << std::endl;
				return false;
		}
		// Check particle velocity
                auto type = (particles.at(ip))->get_type_id();
                vel = (particles.at(ip))->get_velocity();
                if (type == type_ids.at(0)) {
                        if (!(float_equality(vel.at(0), v_mags.at(0)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type A does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(0)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type A does not match expected" << std::endl;
                                return false;
                        }
                } else if (type == type_ids.at(1)) {
                        if (!(float_equality(vel.at(0), v_mags.at(1)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type B does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(1)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type B does not match expected" << std::endl;
                                return false;
                        }

                } else {
                        std::cerr << "Wrong particle type" << std::endl;
                }
	}

	// Check if the values follow uniform distribution
	const double ptcl_x_pos = std::accumulate(all_x_coords.begin(), all_x_coords.end(), 0.0)/Np;
	const double ptcl_y_pos = std::accumulate(all_y_coords.begin(), all_y_coords.end(), 0.0)/Np;
	const double ptcl_heading = std::accumulate(all_headings.begin(), all_headings.end(), 0.0)/Np;
	tol = 0.1;
	if (!(float_equality(ptcl_x_pos, av_x_pos, tol))) {
		std::cerr << "Average x coordinate " << ptcl_x_pos  
				  << " does not match expected " <<  av_x_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_y_pos, av_y_pos, tol))) {
		std::cerr << "Average y coordinate " << ptcl_y_pos  
				  << " does not match expected " <<  av_y_pos << std::endl;
		return false;
	}
	if (!(float_equality(ptcl_heading, av_heading, tol))) {
		std::cerr << "Average heading " << ptcl_heading  
				  << " does not match expected " <<  av_heading << std::endl;
		return false;
	}

	return true;	
	 
}

bool re_initialization_tests_two_types()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	const double tol = 0.1; 	

	// Create the geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	
	// Initialization settings
	const std::vector<size_t> Np_types = {1000, 800};
	const std::vector<double> v_mags = {0.1, 1.0};
	const std::vector<double> R_mins = {0.05, 0.2};
	const std::vector<size_t> type_ids = {0, 1};
	const size_t Np = std::accumulate(Np_types.begin(), Np_types.end(), 0);

	size_t countA = 0;
	size_t countB = 0;

	// Reusable heading
	double theta = 0.0;

	// Initialize and re-initialize particles
	eflow.initialize_particles_from_limits(Np_types, x_lim, y_lim, v_mags, type_ids, R_mins);
	eflow.reinitialize_particles(Np_types, x_lim, y_lim, v_mags, type_ids, R_mins);

	// Collect particles
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();	
	
	// Check particle type count    
        for (size_t ip = 0; ip < particles.size(); ++ ip) {
                if (particles.at(ip)->is_type(0)){
                        ++countA;
                } else if (particles.at(ip)->is_type(1)){
                        ++countB;
                } else {
                        std::cerr << "Error wrong particle type found" << std::endl;
                        return false;
                }
        }

        if (!(countA == Np_types.at(0) && countB == Np_types.at(1))) {
                std::cerr << "Error with assigning particle type" << std::endl;
                return false;
        }


        // Collect positions and headings
        std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
        for (size_t ip = 0; ip < Np; ++ip) {
                coords = (particles.at(ip))->get_position();
                all_x_coords.push_back(coords.at(0));
                all_y_coords.push_back(coords.at(1));
                theta = (particles.at(ip))->get_heading();
                all_headings.push_back(theta);
                // Check limits
                if ((coords.at(0) < x_lim.at(0)) ||
                        (coords.at(0) > x_lim.at(1)) ||
                        (coords.at(1) < y_lim.at(0)) ||
                        (coords.at(1) > y_lim.at(1))) {
                                std::cerr << "x or y coordinates of particle "
                                                  << ip << " are out of bounds" << std::endl;
                                return false;
                }
                // Check particle velocity
                auto type = (particles.at(ip))->get_type_id();
                vel = (particles.at(ip))->get_velocity();
                if (type == type_ids.at(0)) {
                        if (!(float_equality(vel.at(0), v_mags.at(0)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type A does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(0)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type A does not match expected" << std::endl;
                                return false;
                        }
                } else if (type == type_ids.at(1)) {
                        if (!(float_equality(vel.at(0), v_mags.at(1)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type B does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(1)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type B does not match expected" << std::endl;
                                return false;
                        }

                } else {
                        std::cerr << "Wrong particle type" << std::endl;
                }
	}
	
	return true;	
}

bool re_initialization_tests_round_two_types()
{
	// API
	CTC_sim eflow;

	// Floating-point equality threshold
	const double tol = 0.1; 	

	// Create the geometry
	const double R_enclosure = 2.1;
        std::vector<double> enclosure_center = {0.5, 0.1};
        //const double av_x_pos = 0.5, av_y_pos = 0.1, av_heading = 0.0;

	// Initialization settings
	const std::vector<size_t> Np_types = {1000, 800};
	const std::vector<double> v_mags = {0.1, 1.0};
	const std::vector<double> R_mins = {0.05, 0.2};
	const std::vector<size_t> type_ids = {0, 1};
	const size_t Np = std::accumulate(Np_types.begin(), Np_types.end(), 0);

	size_t countA = 0;
	size_t countB = 0;

	// Reusable particle information
        double x0 = 0.0, y0 = 0.0, theta = 0.0;

	// Initialize and re-initialize particles
	eflow.initialize_particles_round_enclosure(Np_types, enclosure_center, R_enclosure, v_mags, type_ids, R_mins);	
	eflow.reinitialize_particles(Np_types, enclosure_center, R_enclosure, v_mags, type_ids, R_mins);


	// Collect particles
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();	
	// Check particle type count    
        for (size_t ip = 0; ip < particles.size(); ++ ip) {
                if (particles.at(ip)->is_type(0)){
                        ++countA;
                } else if (particles.at(ip)->is_type(1)){
                        ++countB;
                } else {
                        std::cerr << "Error wrong particle type found" << std::endl;
                        return false;
                }
        }

        if (!(countA == Np_types.at(0) && countB == Np_types.at(1))) {
                std::cerr << "Error with assigning particle type" << std::endl;
                return false;
        }


       // Check and collect positions and headings
        std::vector<double> all_x_coords, all_y_coords, all_headings, coords, vel;
        for (size_t ip = 0; ip < Np; ++ip) {
                coords = (particles.at(ip))->get_position();
                x0 = coords.at(0);
                y0 = coords.at(1);
                all_x_coords.push_back(x0);
                all_y_coords.push_back(y0);
                theta = (particles.at(ip))->get_heading();
                all_headings.push_back(theta);
                // Check limits
                x0 -= enclosure_center.at(0);
                y0 -= enclosure_center.at(1);
                if ((x0*x0 + y0*y0) > R_enclosure*R_enclosure) {
                                std::cerr << "x or y coordinates of particle "
                                                  << ip << " are out of bounds" << std::endl;
                                return false;
                }
                // Check particle velocity
                auto type = (particles.at(ip))->get_type_id();
                vel = (particles.at(ip))->get_velocity();
                if (type == type_ids.at(0)) {
                        if (!(float_equality(vel.at(0), v_mags.at(0)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type A does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(0)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type A does not match expected" << std::endl;
                                return false;
                        }
                } else if (type == type_ids.at(1)) {
                        if (!(float_equality(vel.at(0), v_mags.at(1)*std::cos(theta), tol))) {
                                std::cerr << "x velocity component of type B does not match expected" << std::endl;
                                return false;
                        }
                        if (!(float_equality(vel.at(1), v_mags.at(1)*std::sin(theta), tol))) {
                                std::cerr << "y velocity component of type B does not match expected" << std::endl;
                                return false;
                        }

                } else {
                        std::cerr << "Wrong particle type" << std::endl;
                }
        }
	
	return true;	
}


// Check application of available boundary conditions
bool boundary_conditions_tests()
{
	if (!make_and_test_boundary_conditions("periodic", "periodic")) 
	{
		std::cerr << "Error with both boundaries periodic" << std::endl;
		return false; 
	}
	if (!make_and_test_boundary_conditions("walls", "walls")) 
	{
		std::cerr << "Error with both boundaries walls" << std::endl;
		return false; 
	}
	if (!make_and_test_boundary_conditions("periodic", "walls")) 
	{
		std::cerr << "Error with x boundaries periodic and y boundaries walls" << std::endl;
		return false; 
	}
	if (!make_and_test_boundary_conditions("walls", "periodic")) 
	{
		std::cerr << "Error with x boundaries walls and y boundaries periodic" << std::endl;
		return false; 
	}
	if (!make_and_test_boundary_conditions_lj("lj-walls", "lj-walls"))
    {
            std::cerr << "Error with both boundaries lj-walls" << std::endl;
            return false;
    }
	if (!make_and_test_boundary_conditions_lj("periodic", "lj-walls"))
    {
            std::cerr << "Error with x boundaries periodic and y boundaries lj-walls" << std::endl;
            return false;
    }
	if (!make_and_test_boundary_conditions_lj("lj-walls", "periodic"))
    {
            std::cerr << "Error with x boundaries lj-walls and y boundaries periodic" << std::endl;
            return false;
    }
	if (!make_and_test_boundary_conditions_lj("walls", "lj-walls"))
    {
            std::cerr << "Error with x boundaries walls and y boundaries lj-walls" << std::endl;
            return false;
    }
    if (!make_and_test_boundary_conditions_lj("lj-walls", "walls"))
    {
            std::cerr << "Error with x boundaries lj-walls and y boundaries walls" << std::endl;
            return false;
    }
	if (!make_and_test_xy_boundary_conditions()) 
	{
		std::cerr << "Error with circular enclosure boundaries" << std::endl;
		return false; 
	}
	return true;	
}

// Setup and test different booundary conditions
bool make_and_test_boundary_conditions_lj(const std::string& x_type, const std::string& y_type)
{
	// API
	CTC_sim eflow;
	

	// Create the geometry (xlim : [0, 1], ylim: [1, 3])
	const double Lx = 1e-6;
	const double Ly = 1e-6;
	const double x_c = 5e-7, y_c = 5e-7;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	const size_t Np = 5;
	const double v_mag = 0.1;
        const double R_particle_wall = 0.25e-6;
	const double epsilon = 0.8, alpha = 0.1, sigma = 1e-7, Fc = -1e-5, dt = 0.1, m_cell = 0.1, tol=1e-5;
	const int m = 9, n = 3; 
	
	// Initialize 5 particles - some outside the domain
	std::vector<double> x_lim = {0.0, Lx};
	std::vector<double> y_lim = {0.0, Ly};
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag, m_cell);
	eflow.set_time_and_noise(dt, 0.0);
	
	// Initial postions (changed by BCs)
	std::vector<std::vector<double>> initial_positions = {
    		{1.01000000e-07, 1.20000000e-07},
    		{2.99000000e-07, 5.00000000e-08},
    		{5.00990000e-07, 1.80000000e-07},
    		{8.99000000e-07, 2.00000000e-08},
    		{7.01000000e-07, 1.00000000e-07}
	};
	// Final positions
	const std::vector<std::vector<double>> exp_lj_lj = {
	    	{4.21208067e+05, 5.83299001e+02},
	    	{-2.89080662e+03, 7.33440029e+08},
	    	{6.07975985e+00, -2.07927432e+04},
	    	{-4.21208067e+05, 7.03110000e+12},
	    	{2.89080662e+03, 4.80036580e+05},
	};

	const std::vector<std::vector<double>> exp_lj_periodic = {
	    	{4.21208067e+05, 1.20000000e-07},
	    	{-2.89080662e+03, 5.00000000e-08},
	    	{6.07975985e+00, 1.80000000e-07},
	    	{-4.21208067e+05, 2.00000000e-08},
		{2.89080662e+03, 1.00000000e-07},
	};

	const std::vector<std::vector<double>> exp_periodic_lj = {
		{1.01000000e-07, 5.83299001e+02},
		{2.99000000e-07, 7.33440029e+08},
		{5.00990000e-07, -2.07927432e+04},
		{8.99000000e-07, 7.03110000e+12},
		{7.01000000e-07, 4.80036580e+05},
	};
	
	// walls–lj  (x clamped to [R_min, L-R_min], y = LJ step)
	const std::vector<std::vector<double>> exp_walls_lj = 
	{
	    	{2.50000000e-07, 5.83299001e+02},
	    	{2.99000000e-07, 7.33440029e+08},
	    	{5.00990000e-07, -2.07927432e+04},
	    	{7.50000000e-07, 7.03110000e+12},
	    	{7.01000000e-07, 4.80036580e+05},
	};

	// lj–walls  (x = LJ step, y clamped to [R_min, L-R_min])
	const std::vector<std::vector<double>> exp_lj_walls = 
	{
	    	{4.21208067e+05, 2.50000000e-07},
	    	{-2.89080662e+03, 2.50000000e-07},
	    	{6.07975985e+00, 2.50000000e-07},
	    	{-4.21208067e+05, 2.50000000e-07},
	    	{2.89080662e+03, 2.50000000e-07},
	};

	// Override each particle position to the known values
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t i = 0; i < Np; ++i) {
        	particles.at(i)->set_new_position(std::move(initial_positions.at(i)));
    }

	// Boundaries 
	if ((x_type == "lj-walls") && (y_type == "lj-walls"))
	{
		// Setup the LJBCs
		eflow.set_all_boundaries_lj(sigma, epsilon, alpha, Fc, m, n, m_cell);	
		// Apply the LJBCs
		eflow.apply_wall_boundary_conditions();

     } else {
		eflow.set_mixed_boundary_conditions(x_type, y_type, R_particle_wall, sigma, epsilon, alpha, Fc, m, n, m_cell);
        eflow.apply_mixed_boundary_conditions(); 
	 }

	// Select correct expected
	const std::vector<std::vector<double>>* expected = nullptr;
	
	if ((x_type == "lj-walls") && (y_type == "lj-walls")) { 
		expected = &exp_lj_lj;
	} else if ((x_type == "lj-walls") && (y_type == "periodic")) {
	    expected = &exp_lj_periodic;
	}
	else if ((x_type == "periodic") && (y_type == "lj-walls")) {
		expected = &exp_periodic_lj;
	}
	else if ((x_type == "walls") && (y_type == "lj-walls")) {
		expected = &exp_walls_lj;
	}
	else if ((x_type == "lj-walls") && (y_type == "walls")) {
        expected = &exp_lj_walls;
	}
	else {
		std::cerr << "Invalid boundary combination: " << x_type << ", " << y_type << std::endl;
    		return false;
	}
	
	// Check if results match
	for (size_t ip = 0; ip < Np; ++ip) {
		auto coords = (particles.at(ip))->get_position();
		if (!(float_equality(coords.at(0), expected->at(ip).at(0), tol) && float_equality(coords.at(1), expected->at(ip).at(1), tol))) {
			std::cerr << "Mismatch at particle " << ip << ": got (" << coords.at(0) 
				<< ", " << coords.at(1) << "), expected ("
                      		<< expected->at(ip).at(0) << ", " << expected->at(ip).at(1) << ")\n";
			return false;
        	}
    	}
	
	return true;
}

// Setup and test different booundary conditions (periodic and lj-walls)
bool make_and_test_boundary_conditions(const std::string& x_type, const std::string& y_type)
{
        // API
        CTC_sim eflow;

        // Create the geometry (xlim : [0, 1], ylim: [1, 3])
        const double Lx = 1.0;
        const double Ly = 2.0;
        const double x_c = 0.5, y_c = 2.0;
        eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

        const size_t Np = 100000;
        const double v_mag = 0.1;
        const double R_particle_wall = 0.25;

        // Initialize 1000 particles - some outside the domain
        std::vector<double> x_lim = {-1.0, 4};
        std::vector<double> y_lim = {-2.5, 6};
        std::vector<double> true_x_lim = {0.0, 1};
        std::vector<double> true_y_lim = {1.0, 3.0};
        std::vector<double> coords;
        eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag);

        // Approach to measure distances
        if ((x_type == "periodic") && (y_type == "periodic"))
        {
                // Setup the PBCs
                eflow.set_all_boundaries_periodic();
                // Apply the PBCs
                eflow.apply_periodic_boundary_conditions();
        } else if ((x_type == "walls") && (y_type == "walls"))
        {
                // Setup the wall BCs
                eflow.set_all_boundaries_walls(R_particle_wall);
                // Apply the wall BCs
                eflow.apply_wall_boundary_conditions();
                // Correct the true limits
                true_x_lim.at(0) += R_particle_wall;
                true_x_lim.at(1) -= R_particle_wall;
                true_y_lim.at(0) += R_particle_wall;
                true_y_lim.at(1) -= R_particle_wall;
        } else {
                eflow.set_mixed_boundary_conditions(x_type, y_type, R_particle_wall);
                eflow.apply_mixed_boundary_conditions();
                // Correct the true limits
                if (x_type == "walls") {
                        true_x_lim.at(0) += R_particle_wall;
                        true_x_lim.at(1) -= R_particle_wall;
                } else {
                        true_y_lim.at(0) += R_particle_wall;
                        true_y_lim.at(1) -= R_particle_wall;
                }
        }

        // Check if all particles are within bounds
        const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
        for (size_t ip = 0; ip < Np; ++ip) {
                coords = (particles.at(ip))->get_position();
                // Check limits
                if ((coords.at(0) < true_x_lim.at(0)) ||
                        (coords.at(0) > true_x_lim.at(1)) ||
                        (coords.at(1) < true_y_lim.at(0)) ||
                        (coords.at(1) > true_y_lim.at(1))) {
                                std::cerr << "x or y coordinates of particle "
                                                  << ip << " are out of bounds" << std::endl;
                                return false;
                }
        }
        return true;
}

// Setup and test different xy boundary conditions (ones that take two arguments)
bool make_and_test_xy_boundary_conditions()
{
	// API
	CTC_sim eflow;
	
	// Create the geometry (xlim : [0, 1], ylim: [1, 3])
	const double radius = 2.0, x_c = 0.5, y_c = 2.0;
	eflow.create_circular_enclosure(radius, x_c, y_c);
	
	const size_t Np = 100000;
	const double v_mag = 0.1;
	const double R_particle_wall = 0.25;
	
	// Initialize 1000 particles - some outside the domain
	std::vector<double> coords;
	eflow.initialize_particles_round_enclosure(Np, {x_c, y_c}, 10.0*radius, v_mag);
	
	// Approach to measure distances
	// Setup the boundaries
	eflow.set_all_boundaries_round_wall(R_particle_wall);
	// Apply the boundary condition 
	eflow.apply_round_wall_boundary_conditions();
	
	// Check if all particles are within bounds
	double x0, y0;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < Np; ++ip) {
		coords = (particles.at(ip))->get_position();
		x0 = coords.at(0);
		y0 = coords.at(1);
	    // Check limits
		x0 -= x_c;
		y0 -= y_c;
		if ((x0*x0 + y0*y0) > radius*radius) {
			std::cerr << "x or y coordinates of particle " 
					  << ip << " are out of bounds" << std::endl;
			return false;
		}
	}
	return true;
}

// Check average heading computations
bool neighbor_heading_tests()
{
	if (!make_and_test_heading("periodic", "periodic")) 
	{
		std::cerr << "Error with both distances periodic" << std::endl;
		return false; 
	}
	if (!make_and_test_heading("walls", "walls")) 
	{
		std::cerr << "Error with both distances walls" << std::endl;
		return false; 
	}
	if (!make_and_test_heading("periodic", "walls")) 
	{
		std::cerr << "Error with x distance periodic and y distance walls" << std::endl;
		return false; 
	}
	if (!make_and_test_heading("walls", "periodic")) 
	{
		std::cerr << "Error with x distance walls and y distance periodic" << std::endl;
		return false; 
	}
	
	if (!make_and_test_heading_two_types("periodic", "periodic"))
	{
		std::cerr << "Error with both distances periodic" << std::endl;
		return false;
	}
	if (!make_and_test_heading_two_types("walls", "walls"))
	{
		std::cerr << "Error with both distances walls" << std::endl;
		return false;
	}
	if (!make_and_test_heading_two_types("periodic", "walls"))
	{
		std::cerr << "Error with x distance periodic and y distance walls" << std::endl;
		return false;
	}
	if (!make_and_test_heading_two_types("walls", "periodic"))
	{
		std::cerr << "Error with x distance walls and y distance periodic" << std::endl;
		return false;
	}
	return true;

}

// Test particle propagation with and without noise
bool particle_propagation_tests()
{
	// API
	CTC_sim eflow;
	
	// Geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	// Particles	
	const size_t Np = 10;
	const double v_mag = 0.1;
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag);

	// Time step and noise
	const double dt = 2.0;
	double eta = 0.0;
	eflow.set_time_and_noise(dt, eta);

	// Test specific
	const size_t max_steps = 2;
	const double tol = 1e-5;
	double new_pos_x = 0.0, new_pos_y = 0.0;
	std::vector<double> old_pos_x(Np), old_pos_y(Np);
    	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();

	// Zero noise test
	std::vector<std::vector<std::vector<double>>> particle_data;
	for (size_t it = 0; it < max_steps; ++it) {
		std::vector<std::vector<double>> current_step_data;
		std::vector<double> headings(Np);
		for (size_t ip = 0; ip < Np; ++ip) {
			old_pos_x.at(ip) = (particles.at(ip))->get_x_position();
			old_pos_y.at(ip) = (particles.at(ip))->get_y_position();
		}
		eflow.propagate_particles();
		for (size_t ip = 0; ip < Np; ++ip) {
			new_pos_x = old_pos_x.at(ip) + v_mag*std::cos((particles.at(ip))->get_heading())*dt;
			new_pos_y = old_pos_y.at(ip) + v_mag*std::sin((particles.at(ip))->get_heading())*dt;
			headings.at(ip) =(particles.at(ip))->get_heading();

			if (!(float_equality(new_pos_x, (particles.at(ip))->get_x_position(), tol))) {
				std::cerr << "Position of particle in x " << (particles.at(ip))->get_x_position()
				  		  << " does not match expected " <<  new_pos_x << std::endl;
				return false;
			}
			if (!(float_equality(new_pos_y, (particles.at(ip))->get_y_position(), tol))) {
				std::cerr << "Position of particle in y " << (particles.at(ip))->get_y_position()
				  		  << " does not match expected " <<  new_pos_y << std::endl;
				return false;
			}
			current_step_data.push_back({new_pos_x, new_pos_y, headings.at(ip)});

		}
		particle_data.push_back(current_step_data);
	}

	print_to_csv(particle_data, "particle_data_zero.csv");

	//  Non-zero noise test
	eta = 0.25;
	eflow.set_time_and_noise(dt, eta);
	std::vector<std::vector<std::vector<double>>> particle_data_non;
	for (size_t it = 0; it < max_steps; ++it) {
		std::vector<std::vector<double>> current_step_data;
		std::vector<double> headings(Np);
		for (size_t ip = 0; ip < Np; ++ip) {
			old_pos_x.at(ip) = (particles.at(ip))->get_x_position();
			old_pos_y.at(ip) = (particles.at(ip))->get_y_position();
		}
		eflow.propagate_particles();
		for (size_t ip = 0; ip < Np; ++ip) {
			new_pos_x = old_pos_x.at(ip) + v_mag*std::cos((particles.at(ip))->get_heading())*dt;
			new_pos_y = old_pos_y.at(ip) + v_mag*std::sin((particles.at(ip))->get_heading())*dt;
			headings.at(ip) =(particles.at(ip))->get_heading();
			if (!(float_equality(new_pos_x, (particles.at(ip))->get_x_position(), tol))) {
				std::cerr << "Position of particle in x " << (particles.at(ip))->get_x_position()
				  		  << " does not match expected " <<  new_pos_x << std::endl;
				return false;
			}
			if (!(float_equality(new_pos_y, (particles.at(ip))->get_y_position(), tol))) {
				std::cerr << "Position of particle in y " << (particles.at(ip))->get_y_position()
				  		  << " does not match expected " <<  new_pos_y << std::endl;
				return false;
			}
			current_step_data.push_back({new_pos_x, new_pos_y, headings.at(ip)});
		}
		particle_data_non.push_back(current_step_data);
	}	

	print_to_csv(particle_data, "particle_data_nonzero.csv");

	return true;
}

/// Setup and test particle interactions 
bool particle_interactions_tests()
{
	if (!make_and_test_interactions("periodic", "periodic")) 
	{
		std::cerr << "Error with both distances periodic" << std::endl;
		return false; 
	}
	if (!make_and_test_interactions("walls", "walls")) 
	{
		std::cerr << "Error with both distances walls" << std::endl;
		return false; 
	}
	if (!make_and_test_interactions("periodic", "walls")) 
	{
		std::cerr << "Error with x distance periodic and y distance walls" << std::endl;
		return false; 
	}
	if (!make_and_test_interactions("walls", "periodic")) 
	{
		std::cerr << "Error with x distance walls and y distance periodic" << std::endl;
		return false; 
	}
	if (!make_and_test_interactions_two_types("periodic", "periodic"))
	{
		std::cerr << "Error with both distances periodic" << std::endl;
		return false;
	}
	if (!make_and_test_interactions_two_types("walls", "walls"))
	{
		std::cerr << "Error with both distances walls" << std::endl;
		return false;
	}
	if (!make_and_test_interactions_two_types("periodic", "walls"))
	{
		std::cerr << "Error with x distance periodic and y distance walls" << std::endl;
		return false;
	}
	if (!make_and_test_interactions_two_types("walls", "periodic"))
	{
		std::cerr << "Error with x distance walls and y distance periodic" << std::endl;
		return false;
	}
	return true;
}

// Sets up and tests interactions computations for different boundary types 
bool make_and_test_interactions(const std::string& x_type, const std::string& y_type)
{
	// API
	CTC_sim eflow;

	// Create the geometry (xlim : [0, 1], ylim: [1, 3])
	const double Lx = 1.0;
	const double Ly = 2.0;
	const double x_c = 0.5, y_c = 2.0;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	const size_t Np = 10;
	const double v_mag = 0.1;
	const double particle_d_min = 0.1;

	// Approach to measure distances 
	if ((x_type == "periodic") && (y_type == "periodic"))
	{
		eflow.set_all_distances_periodic();
	} else if ((x_type == "walls") && (y_type == "walls"))
	{
		eflow.set_all_distances_walls();
	} else {
		eflow.set_mixed_distances(x_type, y_type);
	}

	// Initialize particles
	std::vector<double> x_lim = {0.0, 0.1};
	std::vector<double> y_lim = {1.5, 1.6};	
	std::vector<double> coords;
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag);

	// Initial overlaps
	const int initial_overlaps = count_overlaps(eflow, particle_d_min); 
	if (initial_overlaps == 0) 
	{
		std::cerr << "No overlapping particles" << std::endl;
		return false;
	}

	// Apply particle interactions
	eflow.remove_particle_overlaps(particle_d_min); 
	const int final_overlaps = count_overlaps(eflow, particle_d_min); 

	// Final overlaps (likely 0)
	if (final_overlaps >= initial_overlaps) 
	{
		std::cerr << "Final number of particle overlaps: " << final_overlaps
		 		  << " not smaller than initial: " << initial_overlaps << std::endl;
		return false;
	}
	return true;
}

// Sets up and tests interactions computations for different boundary types 
bool make_and_test_interactions_two_types(const std::string& x_type, const std::string& y_type)
{
	// API
	CTC_sim eflow;

	// Create the geometry (xlim : [0, 1], ylim: [1, 3])
	const double Lx = 1.0;
	const double Ly = 2.0;
	const double x_c = 0.5, y_c = 2.0;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	const std::vector<size_t> Np_types{10, 5};
	const std::vector<double> v_mags = {0.1, 1.0};
        const std::vector<double> R_mins = {0.1, 0.05, 0.2};
        const std::vector<size_t> type_ids = {0, 1};
	

	// Approach to measure distances 
	if ((x_type == "periodic") && (y_type == "periodic"))
	{
		eflow.set_all_distances_periodic();
	} else if ((x_type == "walls") && (y_type == "walls"))
	{
		eflow.set_all_distances_walls();
	} else {
		eflow.set_mixed_distances(x_type, y_type);
	}

	// Initialize particles
	std::vector<double> x_lim = {0.0, 0.1};
	std::vector<double> y_lim = {1.5, 1.6};	
	std::vector<double> coords;
        eflow.initialize_particles_from_limits(Np_types, x_lim, y_lim, v_mags, type_ids, R_mins);
	
	// Initial overlaps
	const int initial_overlaps = count_overlaps(eflow, R_mins);

	if (initial_overlaps == 0) 
	{
		std::cerr << "No overlapping particles" << std::endl;
		return false;
	}

	// Apply particle interactions
	eflow.remove_particle_overlaps_two_types(R_mins); 
	const int final_overlaps = count_overlaps(eflow, R_mins); 

	// Final overlaps (likely 0)
	if (final_overlaps >= initial_overlaps) 
	{
		std::cerr << "Final number of particle overlaps: " << final_overlaps
		 		  << " not smaller than initial: " << initial_overlaps << std::endl;
		return false;
	}
	return true;
}

// Returns the number of overlapping particle pairs
int count_overlaps(CTC_sim& eflow, const double d_min)
{
	int overlaps = 0;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	const size_t Np = particles.size();
	for (size_t ip = 0; ip < Np; ++ip) 
	{
		for (size_t jp = 0; jp < Np; ++jp) 
		{
			if (jp <= ip) 
			{
				continue;
			}
			if (compute_particle_distance(*(particles.at(ip)), *(particles.at(jp))) < d_min)
			{
				++overlaps;
			}
		}
	}
	return overlaps;
}

int count_overlaps(CTC_sim& eflow, const std::vector<double> d_mins)
{
	int overlaps = 0;
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	const size_t Np = particles.size();
	for (size_t ip = 0; ip < Np; ++ip) 
	{
		for (size_t jp = 0; jp < Np; ++jp) 
		{
			if (jp <= ip) 
			{
				continue;
			}
			size_t ti = (particles.at(ip))->get_type_id();
                        size_t tj = (particles.at(jp))->get_type_id();
                        double d_min = d_mins.at(2);
			
			if (ti == tj && ti == 0) {
				d_min = d_mins.at(0);
				if (compute_particle_distance(*(particles.at(ip)), *(particles.at(jp))) < d_min)
				{
					++overlaps;
				}
			} else if (ti == tj && tj == 1) {
				d_min = d_mins.at(1);
                                if (compute_particle_distance(*(particles.at(ip)), *(particles.at(jp))) < d_min)
                                {
                                        ++overlaps;          
                                }	
			} else {
                                if (compute_particle_distance(*(particles.at(ip)), *(particles.at(jp))) < d_min)
                                {
                                        ++overlaps;
                                }
			}
		}
	}
	return overlaps;
}

// Compute distance between two particles
double compute_particle_distance(const Particle& pt1, const Particle& pt2)
{
	const std::vector<double> pos_1 = pt1.get_position();
	const std::vector<double> pos_2 = pt2.get_position();

	return std::sqrt((pos_1.at(0) - pos_2.at(0))*(pos_1.at(0) - pos_2.at(0)) 
		+ (pos_1.at(1) - pos_2.at(1))*(pos_1.at(1) - pos_2.at(1)));
}

// Sets up and tests heading computations for different boundary types
bool make_and_test_heading(const std::string& x_type, const std::string& y_type)
{
	// API
	CTC_sim eflow;
	
	// Geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	// Particles	
	const size_t Np = 100;
	const double v_mag = 0.1;
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v_mag);

	// Approach to measure distances 
	if ((x_type == "periodic") && (y_type == "periodic"))
	{
		eflow.set_all_distances_periodic();
	} else if ((x_type == "walls") && (y_type == "walls"))
	{
		eflow.set_all_distances_walls();
	} else {
		eflow.set_mixed_distances(x_type, y_type);
	} 

	// Averaging
	const size_t max_steps = 2000;
	const double R = 1.5;
	const double tol = 1e-5;
	std::vector<double> old_headings(Np);
	std::vector<double> new_headings(Np);
	std::vector<double> print_headings(Np);

	for (size_t it = 0; it < max_steps; ++it) {
		const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
		for (size_t ip = 0; ip < Np; ++ip) {
			old_headings.at(ip) = (particles.at(ip))->get_heading();
		}
		eflow.compute_average_heading(R);
		for (size_t ip = 0; ip < Np; ++ip) {
			if (ip <= 10) {
				print_headings.at(ip) =(particles.at(ip))->get_heading();
			}
			new_headings.at(ip) = (particles.at(ip))->get_heading();
		}
		if ((it < 30) && (is_equal_floats<double>({old_headings}, {new_headings}, tol))) {
			std::cerr << "Old and new headings should not be equal at step "
					  << it << std::endl;
			return false; 
		}
	
	}

	if (!is_equal_floats<double>({old_headings}, {new_headings}, tol)) {
		std::cerr << "Old and new headings should be equal after "
				  << max_steps << " steps " << std::endl;
		return false; 
	}

	std::ofstream outputFile("particle_headings_" + x_type + "_" + y_type + ".csv");
	for (size_t i = 0; i < print_headings.size(); ++i) 
	{
		outputFile << print_headings.at(i);
		if (i < print_headings.size() - 1) {
			outputFile << ",";
		}
	}
	outputFile.close();

	return true;
}

// Sets up and tests heading computations for different boundary types (two particle types)
bool make_and_test_heading_two_types(const std::string& x_type, const std::string& y_type)
{
	// API
	CTC_sim eflow;

	// Geometry
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	
	// Particles
	const std::vector<size_t> Np_types = {1000, 800};
        const std::vector<double> v_mags = {0.1, 1.0};
        const std::vector<double> R_mins = {0.05, 0.2};
        const std::vector<size_t> type_ids = {0, 1};
        const size_t Np = std::accumulate(Np_types.begin(), Np_types.end(), 0);
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
        eflow.initialize_particles_from_limits(Np_types, x_lim, y_lim, v_mags, type_ids, R_mins);

	// Approach to measure distances
	if ((x_type == "periodic") && (y_type == "periodic"))
	{
		eflow.set_all_distances_periodic();
	} else if ((x_type == "walls") && (y_type == "walls"))
	{
		eflow.set_all_distances_walls();
	} else {
		eflow.set_mixed_distances(x_type, y_type);
	}

	// Averaging
	const size_t max_steps = 2000;
	const std::vector<double> R = {1.5, 3.0};
	const double tol = 1e-5;
	std::vector<double> old_headings(Np);
	std::vector<double> new_headings(Np);
	std::vector<double> print_headings(Np);

	for (size_t it = 0; it < max_steps; ++it) {
		const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
		for (size_t ip = 0; ip < Np; ++ip) {
			old_headings.at(ip) = (particles.at(ip))->get_heading();
		}
		eflow.compute_average_heading_two_types(R, false);
		for (size_t ip = 0; ip < Np; ++ip) {
			if (ip <= 10) {
				print_headings.at(ip) =(particles.at(ip))->get_heading();
			}
			new_headings.at(ip) = (particles.at(ip))->get_heading();
		}
		if ((it < 30) && (is_equal_floats<double>({old_headings}, {new_headings}, tol))) {
			std::cerr << "Old and new headings should not be equal at step "
					  << it << std::endl;
			return false;
		}

	}

	if (!is_equal_floats<double>({old_headings}, {new_headings}, tol)) {
		std::cerr << "Old and new headings should be equal after "
				  << max_steps << " steps " << std::endl;
		return false;
	}

	std::ofstream outputFile("mixed_particle_headings_" + x_type + "_" + y_type + ".csv");
	for (size_t i = 0; i < print_headings.size(); ++i)
	{
		outputFile << print_headings.at(i);
		if (i < print_headings.size() - 1) {
			outputFile << ",";
		}
	}
	outputFile.close();

	return true;
}

// Check if the program reacts to wrong types of boundary conditions
bool test_wrong_boundary_types()
{
	// API
	CTC_sim eflow;
	// For wall BCs
	const double R_particle_wall = 0.25;
	// For lj-walls BCs
	const double epsilon = 0.8, alpha = 0.1, sigma = 1e-7, Fc = -1e-5, m_cell = 0.1;
        const int m = 9, n = 3;


	// Settings for tests
	bool verbose = true;
	const std::string x_type{"periodic"}, x_wrong_type{"pressure"}, y_type{"walls"}, y_wrong_type{"outlet"}, lj_x_wrong_type{"ljwalls"}, lj_y_wrong_type{"ljwalls"};
	const std::invalid_argument invarg_x("Wrong type of x-boundary: " + x_wrong_type);	
	const std::invalid_argument invarg_y("Wrong type of y-boundary: " + y_wrong_type);

	// Wrong in x, then wrong in y
	/*if (!exception_test(verbose, &invarg_x, &CTC_sim::set_mixed_boundary_conditions, eflow, x_wrong_type, y_type, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, R_particle_wall))
	{
		std::cerr << "Wrong x boundary type not recognized as an error" << std::endl;
		return false;
	}
	if (!exception_test(verbose, &invarg_y, &CTC_sim::set_mixed_boundary_conditions, eflow, x_type, y_wrong_type, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, R_particle_wall))
	{
		std::cerr << "Wrong y boundary type not recognized as an error" << std::endl;
		return false;
	}
	// Test for LJ
	if (!exception_test(verbose, &invarg_x, &CTC_sim::set_mixed_boundary_conditions, eflow, lj_x_wrong_type, y_type, sigma, epsilon, alpha, Fc, m, n, m_cell, 0.0))
        {
                std::cerr << "Wrong x boundary type (LJ typo) not recognized as an error" << std::endl;
                return false;
        }
        if (!exception_test(verbose, &invarg_y, &CTC_sim::set_mixed_boundary_conditions, eflow, x_wrong_type, y_type, sigma, epsilon, alpha, Fc, m, n, m_cell, 0.0))
        {
                std::cerr << "Wrong y boundary type (LJ typo) not recognized as an error" << std::endl;
                return false;
        }
*/

	return true;
}

// Check if the program reacts to wrong types of distance computation approaches 
bool test_wrong_distance_types()
{
	// API
	CTC_sim eflow;

	// Settings for tests
	bool verbose = true;
	const std::string x_type{"periodic"}, x_wrong_type{"zero"}, y_type{"walls"}, y_wrong_type{"gradient"};
	const std::invalid_argument invarg_x("Wrong type of x-boundary: " + x_wrong_type);	
	const std::invalid_argument invarg_y("Wrong type of y-boundary: " + y_wrong_type);

	// Wrong in x, then wrong in y
	if (!exception_test(verbose, &invarg_x, &CTC_sim::set_mixed_distances, eflow, x_wrong_type, y_type))
	{
		std::cerr << "Wrong x distance computation type not recognized as an error" << std::endl;
		return false;
	}
	if (!exception_test(verbose, &invarg_y, &CTC_sim::set_mixed_distances, eflow, x_type, y_wrong_type))
	{
		std::cerr << "Wrong y distance computation type not recognized as an error" << std::endl;
		return false;
	}
	return true;
}

// Test suite for introduction of the external flow
bool external_flow()
{
	// Results and comparison
	std::vector<double> temp_data;
	const double tol = 1e-5;
	// Flow field
	// Needs to match the value in apply_external_flow
	const double Ly = 2.0;
	// External velocity magnitude
	double u_mag = 0.0;
	// Vicsek velocity magnitude
	double u_vic = 0.0;
	// Everything with "flow_test_" prefix is a static variable
	// Defined at the top level of this file
	flow_test_channel_radius_squared = Ly*Ly/4;
	flow_test_channel_radius = Ly/2;
	flow_test_mu = 1e-3;

	// x component (uses static variables)
	auto vel_x = [](double x, double y) 
		{ return (1.0/(4.0*flow_test_mu)*flow_test_dpdL*
					(flow_test_channel_radius_squared - (y - flow_test_channel_radius)*
					(y - flow_test_channel_radius))); }; 

	// y component
	auto vel_y = [](double x, double y) { return 0.0; };

	// Pure Vicsek
	u_mag = 0.0;
	u_vic = 0.5;
	std::vector<std::vector<double>> pure_vicsek = apply_external_flow(u_vic, 0.0, 0.2);
	// The velocities of particles should be smaller or equal to Vicsek velocity magnitude
	for (size_t jf = 0; jf < pure_vicsek.size(); ++jf) 
	{
		temp_data = pure_vicsek.at(jf);
		if (temp_data.at(2) > u_vic) 
		{
			std::cerr << "x velocity component of particle " << jf << " " << temp_data.at(2)			  		  
				<< " larger than expected " <<  u_vic << std::endl;
			return false;
		}
		if (temp_data.at(3) > u_vic) 
		{
			std::cerr << "y velocity component of particle " << jf << " " << temp_data.at(3)			  		  
				<< " larger than expected " <<  u_vic << std::endl;
			return false;
		}
	}

	// Pure flow (particles as passive tracers with self alignment)
	u_mag = 2.0;
	u_vic = 0.0;
	std::vector<std::vector<double>> pure_flow = apply_external_flow(u_vic, u_mag, 0.0);
	flow_test_dpdL = 8*flow_test_mu*u_mag/flow_test_channel_radius_squared;
	// The velocities of particles should be equal to flow velocities
	for (size_t jf = 0; jf < pure_flow.size(); ++jf) 
	{
		temp_data = pure_flow.at(jf);
		if (!(float_equality(vel_x(temp_data.at(0), temp_data.at(1)), temp_data.at(2), tol)))
		{
			std::cerr << "x velocity component of particle " << jf << " " << temp_data.at(2)			  		  
				<< " does not match the external flow " <<  vel_x(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
		if (!(float_equality(vel_y(temp_data.at(0), temp_data.at(1)), temp_data.at(3), tol)))
		{
			std::cerr << "y velocity component of particle " << jf << " " << temp_data.at(3)			  		  
				<< " does not match the external flow " <<  vel_y(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
	}

	// Both flows - particle velocity much smaller than external
	u_mag = 2.0;
	u_vic = 0.1;
	std::vector<std::vector<double>> vicsek_smaller = apply_external_flow(u_vic, u_mag, 0.2);
	flow_test_dpdL = 8*flow_test_mu*u_mag/flow_test_channel_radius_squared;
	// The velocities of particles should not be equal to flow velocities
	// with an exception of zero
	for (size_t jf = 0; jf < pure_flow.size(); ++jf) 
	{
		temp_data = vicsek_smaller.at(jf);
		if ((float_equality(vel_x(temp_data.at(0), temp_data.at(1)), temp_data.at(2), tol))
			&& !(float_equality(0.0, temp_data.at(2), tol))) 
		{
			std::cerr << "x velocity component of particle " << jf << " " << temp_data.at(2)			  		  
				<< " matches the external flow " <<  vel_x(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
		if ((float_equality(vel_y(temp_data.at(0), temp_data.at(1)), temp_data.at(3), tol))
	 		&& !(float_equality(0.0, temp_data.at(3), tol)))
		{
			std::cerr << "y velocity component of particle " << jf << " " << temp_data.at(3)			  		  
				<< " matches the external flow " <<  vel_y(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
	}

	// Both flows - particle velocity similar external
	u_mag = 2.0;
	u_vic = 2.0;
	std::vector<std::vector<double>> vicsek_equal = apply_external_flow(u_vic, u_mag, 0.2);
	flow_test_dpdL = 8*flow_test_mu*u_mag/flow_test_channel_radius_squared;
	// The velocities of particles should not be equal to flow velocities
	// with an exception of zero
	for (size_t jf = 0; jf < pure_flow.size(); ++jf) 
	{
		temp_data = vicsek_equal.at(jf);
		if ((float_equality(vel_x(temp_data.at(0), temp_data.at(1)), temp_data.at(2), tol))
			&& !(float_equality(0.0, temp_data.at(2), tol))) 
		{
			std::cerr << "x velocity component of particle " << jf << " " << temp_data.at(2)			  		  
				<< " matches the external flow " <<  vel_x(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
		if ((float_equality(vel_y(temp_data.at(0), temp_data.at(1)), temp_data.at(3), tol))
	 		&& !(float_equality(0.0, temp_data.at(3), tol)))
		{
			std::cerr << "y velocity component of particle " << jf << " " << temp_data.at(3)			  		  
				<< " matches the external flow " <<  vel_y(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
	}

	// Both flows - particle velocity much larger than external
	u_mag = 0.1;
	u_vic = 2.0;
	std::vector<std::vector<double>> vicsek_higher = apply_external_flow(u_vic, u_mag, 0.2);
	flow_test_dpdL = 8*flow_test_mu*u_mag/flow_test_channel_radius_squared;
	// The velocities of particles should not be equal to flow velocities
	// with an exception of zero
	for (size_t jf = 0; jf < pure_flow.size(); ++jf) 
	{
		temp_data = vicsek_higher.at(jf);
		if ((float_equality(vel_x(temp_data.at(0), temp_data.at(1)), temp_data.at(2), tol))
			&& !(float_equality(0.0, temp_data.at(2), tol))) 
		{
			std::cerr << "x velocity component of particle " << jf << " " << temp_data.at(2)			  		  
				<< " matches the external flow " <<  vel_x(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
		if ((float_equality(vel_y(temp_data.at(0), temp_data.at(1)), temp_data.at(3), tol))
	 		&& !(float_equality(0.0, temp_data.at(3), tol)))
		{
			std::cerr << "y velocity component of particle " << jf << " " << temp_data.at(3)			  		  
				<< " matches the external flow " <<  vel_y(temp_data.at(0), temp_data.at(1)) << std::endl;
			return false;
		}
	}

	return true;
}

// Test for external force application
bool external_force()
{
	// --- Compare the computed with calculated positions

	// Overwrite the static variables
	F_x_max = 100.0;		
	F_y_max = 20.0;
	step_size = 3.0;
	domain_height = 10.0;

	// Other variables
	double mass = 10.0;
	double dt = 0.5;
	
	// x component (uses static variables)
    auto F_x = [](double y) 
	{
		if (((y >= 0) && (y <= step_size)) || ((y >= domain_height - step_size) && (y <= domain_height)))
		{
   			return F_x_max;
		} else {
			return 0.0;
		}
    };
	
	// y component (uses static variables)
	auto F_y = [](double y) 
	{
		if (((y >= 0) && (y <= step_size)) || ((y >= domain_height - step_size) && (y <= domain_height)))
		{
   			return F_y_max;
		} else {
			return 0.0;
		}
    };

	auto vel_anything = [](double x, double y) { return 0.0; };
	
	CTC_sim eflow(vel_anything, vel_anything, F_x, F_y);

	// Channel dimensions 
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;
	// Geometry	
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);
	// Particles	
	const size_t Np = 10;
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, 0, mass);
	// Time step and noise
	eflow.set_time_and_noise(dt, 0);

	// Calculate expected positions 
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	std::vector<double> results_x;
	std::vector<double>	results_y;
	for (size_t ip = 0; ip < Np; ++ip) 
	{
		double x = (particles.at(ip))->get_x_position();
		double y = (particles.at(ip))->get_y_position();
		results_x.push_back((x + F_x(y)*dt*dt/mass));
		results_y.push_back((y + F_y(y)*dt*dt/mass));
	}

	// Compute using API
	eflow.apply_external_force_field();

	// Compare
	double tol = 1e-5;
	for (size_t ip = 0; ip < Np; ++ip) 
	{
		double x = (particles.at(ip))->get_x_position();
		double y = (particles.at(ip))->get_y_position();
		double xexp = results_x.at(ip);
		double yexp = results_y.at(ip);
		if (!(float_equality(x, xexp, tol))) {
			std::cerr << "Position of particle in x " << x 
			  		  << " does not match expected " <<  xexp << std::endl;
			return false;
		}
		if (!(float_equality(y, yexp, tol))) {
			std::cerr << "Position of particle in y " << y 
			  		  << " does not match expected " <<  yexp << std::endl;
			return false;
		}
	}

	return true;
}

// Tests bonds, lift force, and new computation (from a single reference point)
bool lift_force_test()
{
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

	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Time stepping and collection settings
	// Time step
	double dt = particle_parameters.at("Time step");
	// Max simulation time
	size_t max_steps = static_cast<size_t>(particle_parameters.at("Maximum simulation time")/dt);

	// Collection times (time 0 is collected by default)
	std::vector<size_t> save_times(2);
	std::iota(save_times.begin(), save_times.end(), 1.0);
	// Convert to steps
	std::vector<size_t> save_steps_vector(save_times.size());
	std::transform(save_times.cbegin(), save_times.cend(), save_steps_vector.begin(), [&dt](size_t ti) 
		{ return static_cast<size_t>(static_cast<double>(ti)/dt); });
	std::set<size_t> save_steps(save_steps_vector.begin(), save_steps_vector.end());

	// Also noise
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));
	// And new variables
	eflow.set_lift_force_coefficient(1.53, 3.8);

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

			// Compute bond forces, check for failure, add contributions
			eflow.process_bonds();			

			// Lift force
			eflow.add_lift_force(particle_parameters.at("Particle radius"),
				particle_parameters.at("Density of the fluid"),
				particle_parameters.at("Viscosity of the fluid"));

			// Save active bond count
			eflow.track_bonds();

			// Compare to bond pairs
			current_bonds_with_time = eflow.get_bond_number_with_time();
			if (current_bonds_with_time.back() != bond_pairs.size())
			{
				std::cerr << "Number of bonds not matching number of ID pairs" << std::endl;
			}
		
			// Update everything for this step
			eflow.update_positions_and_headings();
	
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

	return true;	
}

// Simulate particles with external flow in a 2D channel
// Arguments are particle and external flow velocity magnitudes (v and u) and noise (eta)
// Returns particle coordinates and velocities
std::vector<std::vector<double>> apply_external_flow(const double v, const double u, const double eta)
{
	// Final x and y coordinates and velocity components
	// of all particles
	std::vector<std::vector<double>> results;
	std::vector<double> temp_velocity{0, 0};
	std::vector<double> temp_results{0, 0, 0, 0};

	// Channel dimensions 
	const double Lx = 10.0;
	const double Ly = 2.0;
	const double x_c = 1.0, y_c = 0.1;

	// External flow - laminar flow through a 2D channel
	// Everything with "flow_test_" prefix is a static variable
	// Defined at the top level of this file
	flow_test_channel_radius_squared = Ly*Ly/4;
	flow_test_channel_radius = Ly/2;
	flow_test_mu = 1e-3;
	flow_test_dpdL = 8*flow_test_mu*u/flow_test_channel_radius_squared;

	const bool needs_external_flow = true;

	// x component (uses static variables)
	auto vel_x = [](double x, double y) 
		{ return (1.0/(4.0*flow_test_mu)*flow_test_dpdL*
					(flow_test_channel_radius_squared - (y - flow_test_channel_radius)*
					(y - flow_test_channel_radius))); }; 

	// y component
	auto vel_y = [](double x, double y) { return 0.0; };

	// API
	CTC_sim eflow(vel_x, vel_y);

	// Geometry	
	eflow.create_rectangular_enclosure(Lx, Ly, x_c, y_c);

	// Boundary conditions and distances
	std::string x_type{"periodic"}, y_type{"walls"};
	// Minimum distance from a wall
	const double R_particle_wall = 1e-8;

	eflow.set_mixed_distances(x_type, y_type);
	eflow.set_mixed_boundary_conditions(x_type, y_type, R_particle_wall);
	eflow.apply_mixed_boundary_conditions();

	// Particles	
	const size_t Np = 10;
	// Neighbor definition
	const double R_neighbor = 0.5;
	std::vector<double> all_lim = eflow.get_enclosure_limits();
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};

	eflow.initialize_particles_from_limits(Np, x_lim, y_lim, v);

	// Time step and noise
	const size_t max_steps = 100;
	const double dt = 2.0;

	eflow.set_time_and_noise(dt, eta);

	// Simulation
	for (size_t it = 0; it < max_steps; ++it) 
	{
		eflow.compute_average_heading(R_neighbor, needs_external_flow);
		eflow.propagate_particles(needs_external_flow);
		eflow.apply_mixed_boundary_conditions();
	}

	// Collect results
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();
	for (size_t ip = 0; ip < Np; ++ip) 
	{
		temp_results.at(0) = (particles.at(ip))->get_x_position();
		temp_results.at(1) = (particles.at(ip))->get_y_position();
		temp_velocity = (particles.at(ip))->get_velocity();
		temp_results.at(2) = temp_velocity.at(0);
		temp_results.at(3) = temp_velocity.at(1);
		results.push_back(temp_results);
	}

	return results;
	
}

// Tests the correct application of the lift force
bool lift_force_contribution_test()
{
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

	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Also noise
	double dt = 10;
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));
	// And new variables
	eflow.set_lift_force_coefficient(20.0, 30.0);

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
									- 2.0*particle_parameters.at("Particle radius"),
								particle_parameters.at("Domain center y")
									- 2.0*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center x")
									+ 2*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center y")
									+ 2*particle_parameters.at("Particle radius") };
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.set_initial_heading_interval(-EASY_PI/2.0, EASY_PI/2.0);
	
	// Initialize particles and bonds
	double mass = 7;
	eflow.initialize_particles_from_limits(particle_parameters.at("Number of particles"), 
			x_lim, y_lim, particle_parameters.at("Velocity magnitude")*1e5, mass);
	eflow.install_bonds(particle_parameters.at("Bond normal stiffness"), 
						particle_parameters.at("Bond shear stiffness"),
						particle_parameters.at("Bond normal strength"),
						particle_parameters.at("Bond shear strength"),
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap")/3);

	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);
	// For particle properties
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();

	// Collect the initial x, y coordinates of particles 
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		temp_info.at(0) = particles.at(ip)->get_x_position();
		temp_info.at(1) = particles.at(ip)->get_y_position();
		temp_info.at(2) = (particles.at(ip))->get_x_velocity();
		temp_info.at(3) = static_cast<double>(0);
		particle_info.push_back(temp_info);
	}

	// Lift force
	eflow.add_lift_force(0.1, 1e3, 1e-3);
	eflow.update_positions_and_headings();

	// Compare with expected displacement
	double tol = 1e-5, F_lift = 0.0;
	const double dy_step = 1.0e-6; 
	double pre_cell_F = 1.292, pre_cluster_F = 1.938;
	for (size_t ip = 0; ip < particles.size(); ++ip) 
	{
		double x = particle_info.at(ip).at(0);
		double y = particle_info.at(ip).at(1);
		double u_particle = particle_info.at(ip).at(2);
		
		// Slip velocity
		// 2. Get Fluid Velocity (Query the External Flow Object directly)
        double u_fluid = vel_x(x, y);

        // 3. Calculate Local Shear Rate (du/dy) numerically
        // We sample velocity slightly above and below to find the slope
        double u_upper = vel_x(x, y + dy_step);
        double u_lower = vel_x(x, y - dy_step);
        
        // Central Difference Approximation: slope = (f(x+h) - f(x-h)) / 2h
        double shear_rate = (u_upper - u_lower) / (2.0 * dy_step);
        double abs_shear = std::fabs(shear_rate);

        // 4. Calculate Slip Velocity
        double u_slip = u_fluid - u_particle;

        // 5. Compute Lift Force Magnitude
        // Formula: F ~ 6.46 * R^2 * sqrt(rho * mu * |shear|) * u_slip
		if (static_cast<bool>(particles.at(ip)->get_bond_number()))
		{
        	F_lift = pre_cluster_F * std::sqrt(abs_shear) * u_slip;
			//std::cout << "Some are clusters" << std::endl;
		} else {
			F_lift = pre_cell_F * std::sqrt(abs_shear) * u_slip;
			//std::cout << "Some are single cells" << std::endl;
		}

        // 6. Apply Direction
        // If shear is positive (velocity increases with Y), lagging particles go UP (+Y)
        // If shear is negative (velocity decreases with Y), lagging particles go DOWN (-Y)
        double direction_sign = (shear_rate > 0) ? 1.0 : -1.0;
        double force_y = F_lift * direction_sign;

		// Expected positions
		double xexp = particle_info.at(ip).at(0);
		double yexp = particle_info.at(ip).at(1) + force_y*dt*dt/mass;

		x =  particles.at(ip)->get_x_position();
		y =  particles.at(ip)->get_y_position();

		if (!(float_equality(x, xexp, tol))) {
			std::cerr << "Position of particle in x " << x 
			  		  << " does not match expected " <<  xexp << std::endl;
			return false;
		}
		if (!(float_equality(y, yexp, tol))) {
			std::cerr << "Position of particle in y " << y 
			  		  << " does not match expected " <<  yexp << std::endl;
			return false;
		}
	}
	return true;	
}

// Tests the correct addition of bonds and lift contributions 
bool all_contributions_test()
{
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

	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Also noise
	double dt = 10;
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));
	// And new variables
	eflow.set_lift_force_coefficient(200.0, 300.0);

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
									- 2.0*particle_parameters.at("Particle radius"),
								particle_parameters.at("Domain center y")
									- 5.0*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center x")
									+ 2*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center y")
									+ 2*particle_parameters.at("Particle radius") };
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.set_initial_heading_interval(-EASY_PI/2.0, EASY_PI/2.0);
	
	// Initialize particles and bonds
	double mass = 7;
	eflow.initialize_particles_from_limits(particle_parameters.at("Number of particles"), 
			x_lim, y_lim, particle_parameters.at("Velocity magnitude"), mass);
	eflow.install_bonds(particle_parameters.at("Bond normal stiffness")*1e14, 
						particle_parameters.at("Bond shear stiffness")*1e14,
						particle_parameters.at("Bond normal strength")*1e16,
						particle_parameters.at("Bond shear strength")*1e16,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));
	
	eflow.install_wall_bonds(particle_parameters.at("Wall bond normal stiffness")*1e10, 
						particle_parameters.at("Wall bond shear stiffness")*1e10,
						particle_parameters.at("Wall bond normal strength")*1e10,
						particle_parameters.at("Wall bond shear strength")*1e10,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Domain length y")*0.45);

	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);
	// For particle properties
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();

	// Collect the initial x, y coordinates of particles 
	for (size_t ip = 0; ip < particles.size(); ++ip) {
		temp_info.at(0) = particles.at(ip)->get_x_position();
		temp_info.at(1) = particles.at(ip)->get_y_position();
		temp_info.at(2) = static_cast<double>(0);
		temp_info.at(3) = static_cast<double>(0);
		particle_info.push_back(temp_info);
	}

	for (size_t ir = 0; ir < 2; ++ir)
	{
		// Lift force
		eflow.add_lift_force(0.1, 1e3, 1e-3);
		
		// Align heading with neighboring cells
		eflow.compute_average_heading(particle_parameters.at("Interaction radius"),
				needs_external_flow);
	
		// Update particle positions
		eflow.propagate_particles(needs_external_flow);
	
		// Compute bond forces, check for failure, add contributions
		eflow.process_bonds();		

		// Compute bond forces, check for failure, add contributions
		eflow.process_wall_bonds();		
		
		// Collect contributions
		std::vector<std::vector<double>> all_x;
		std::vector<std::vector<double>> all_y; 
		std::vector<std::vector<double>> delta_headings;
		std::vector<double>	new_headings;

		for (size_t ip = 0; ip < particles.size(); ++ip) 
		{
 			all_x.push_back((particles.at(ip))->get_all_x_contributions()); 
			all_y.push_back((particles.at(ip))->get_all_y_contributions());
			new_headings.push_back((particles.at(ip))->get_new_heading());
			delta_headings.push_back((particles.at(ip))->get_heading_increments());
		}

/*		for (size_t ip = 0; ip < particles.size(); ++ip) 
		{
		    std::cout << "Particle " << ip << ":\n";
		
		    std::cout << "  all_x: ";
		    for (size_t j = 0; j < all_x.at(ip).size(); ++j) 
		    {
		        std::cout << all_x.at(ip).at(j) << " ";
		    }
		    std::cout << "\n";
		
		    std::cout << "  all_y: ";
		    for (size_t j = 0; j < all_y.at(ip).size(); ++j) 
		    {
		        std::cout << all_y.at(ip).at(j) << " ";
		    }
		    std::cout << "\n";
		
		    std::cout << "  delta_headings: ";
		    for (size_t j = 0; j < delta_headings.at(ip).size(); ++j) 
		    {
		        std::cout << delta_headings.at(ip).at(j) << " ";
		    }
		    std::cout << "\n";
		
		    std::cout << "  new_heading: " << new_headings.at(ip) << "\n\n";
		}
*/
		eflow.update_positions_and_headings();
		                                                                   	
		// Compare with expected displacement
		double tol = 1e-4;
		for (size_t ip = 0; ip < particles.size(); ++ip) 
		{
			double x = (particles.at(ip))->get_x_position();
			double y = (particles.at(ip))->get_y_position();
			double theta = (particles.at(ip))->get_heading();
			double xexp = particle_info.at(ip).at(0);
			double yexp = particle_info.at(ip).at(1);
			double theta_exp = new_headings.at(ip) + std::accumulate(delta_headings.at(ip).begin(), delta_headings.at(ip).end(), 0.0);
		
			xexp += std::accumulate(all_x.at(ip).begin(), all_x.at(ip).end(), 0.0);
			yexp += std::accumulate(all_y.at(ip).begin(), all_y.at(ip).end(), 0.0);

			if (!(float_equality(x, xexp, tol))) {
				std::cerr << "Position of particle in x " << x 
				  		  << " does not match expected " <<  xexp << " at step " << ir << std::endl;
				return false;
			}
			if (!(float_equality(y, yexp, tol))) {
				std::cerr << "Position of particle in y " << y 
				  		  << " does not match expected " <<  yexp << " at step " << ir << std::endl;
				return false;
			}
			if (!(float_equality(theta, theta_exp, tol))) {
				std::cerr << "Particle heading " << theta 
				  		  << " does not match expected " <<  theta_exp << " at step " << ir << std::endl;
				return false;
			}
		}

		particle_info.clear();

		// Collect the initial x, y coordinates of particles 
		for (size_t ip = 0; ip < particles.size(); ++ip) {
			temp_info.at(0) = particles.at(ip)->get_x_position();
			temp_info.at(1) = particles.at(ip)->get_y_position();
			temp_info.at(2) = static_cast<double>(0);
			temp_info.at(3) = static_cast<double>(0);
			particle_info.push_back(temp_info);
		}

	}
	return true;	
}

// Tests the new distance-based bond breaking logic 
/*bool dist_bond_breaking_test()
{
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

	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Also noise
	double dt = 10;
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));
	// And new variables
	eflow.set_lift_force_coefficient(200.0, 300.0);

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
									- 2.0*particle_parameters.at("Particle radius"),
								particle_parameters.at("Domain center y")
									- 5.0*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center x")
									+ 2*particle_parameters.at("Particle radius"),
							    particle_parameters.at("Domain center y")
									+ 2*particle_parameters.at("Particle radius") };
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.set_initial_heading_interval(-EASY_PI/2.0, EASY_PI/2.0);
	
	// Initialize particles and bonds
	double mass = 7;
	// Number of particles
	size_t Np = 200;
	eflow.initialize_particles_from_limits(Np, 
			x_lim, y_lim, particle_parameters.at("Velocity magnitude"), mass);
	eflow.install_bonds(particle_parameters.at("Bond normal stiffness")*1e14, 
						particle_parameters.at("Bond shear stiffness")*1e14,
						particle_parameters.at("Bond normal strength")*1e16,
						particle_parameters.at("Bond shear strength")*1e16,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"), particle_parameters.at("Max bond distance"));
	
	eflow.install_wall_bonds(particle_parameters.at("Wall bond normal stiffness")*1e10, 
						particle_parameters.at("Wall bond shear stiffness")*1e10,
						particle_parameters.at("Wall bond normal strength")*1e10,
						particle_parameters.at("Wall bond shear strength")*1e10,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Domain length y")*0.45);

	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);
	// For particle properties
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();

	// Bond info
	std::list<Bond>& bonds = eflow.get_particle_bonds();

	for (size_t ir = 0; ir < 200; ++ir)
	{
		// Lift force
		eflow.add_lift_force(0.1, 1e3, 1e-3);
		
		// Align heading with neighboring cells
		eflow.compute_average_heading(particle_parameters.at("Interaction radius"),
				needs_external_flow);
	
		// Update particle positions
		eflow.propagate_particles(needs_external_flow);
	
		// Compute bond forces, check for failure, add contributions
		eflow.process_bonds();		

		// Compute bond forces, check for failure, add contributions
		eflow.process_wall_bonds();		

		for (auto& bond_i : bonds)
		{
			// Get the IDs of bonded particles
			std::vector<size_t> pIDs = bond_i.get_particle_IDs();
	
			// Compute particle distance
			Particle temp_particle_i = *(particles.at(pIDs.at(0)));
			Particle temp_particle_j = *(particles.at(pIDs.at(1)));
	
			std::vector<double> temp_coords = temp_particle_i.get_position();
			double xi = temp_coords.at(0); 			
			double yi = temp_coords.at(1);
			temp_coords = temp_particle_j.get_position();
			double xj = temp_coords.at(0); 			
			double yj = temp_coords.at(1);
		
			// Distance (x - periodic, y - walls)
			double dx = std::min<double>(std::fabs(xi - xj), 
					particle_parameters.at("Domain length x") - std::fabs(xi - xj));
			double dy = yi - yj;
		
			double dsq = dx*dx + dy*dy;
	
			// Break the bond if distance is too large
			if (std::sqrt(dsq) > particle_parameters.at("Max bond distance"))
			{
				std::cerr << "Long bond not removed " << std::endl;
				return false;			
			}
		}

		eflow.update_positions_and_headings();

		eflow.reinstall_bonds(particle_parameters.at("Bond normal stiffness")*1e14, 
						particle_parameters.at("Bond shear stiffness")*1e14,
						particle_parameters.at("Bond normal strength")*1e16,
						particle_parameters.at("Bond shear strength")*1e16,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));

	}
	return true;	
}

// Tests the new distance-based wall bond breaking logic 
bool dist_wall_bond_breaking_test()
{
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

	// API
	CTC_sim eflow(vel_x, vel_y);
	
	// Also noise
	double dt = 10;
	eflow.set_time_and_noise(dt, particle_parameters.at("Amplitude of noise"));
	// And new variables
	eflow.set_lift_force_coefficient(20.0, 30.0);

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
									- 2.0*particle_parameters.at("Particle radius"),
								0.0,
							    particle_parameters.at("Domain center x")
									+ 2*particle_parameters.at("Particle radius"),
								2*particle_parameters.at("Particle radius") };
	std::vector<double> x_lim{all_lim.at(0), all_lim.at(2)};
	std::vector<double> y_lim{all_lim.at(1), all_lim.at(3)};
	eflow.set_initial_heading_interval(-EASY_PI/2.0, EASY_PI/2.0);
	
	// Initialize particles and bonds
	double mass = 7;
	// Number of particles
	size_t Np = 200;
	eflow.initialize_particles_from_limits(Np, 
			x_lim, y_lim, particle_parameters.at("Velocity magnitude"), mass);
	eflow.install_bonds(particle_parameters.at("Bond normal stiffness")*1e14, 
						particle_parameters.at("Bond shear stiffness")*1e14,
						particle_parameters.at("Bond normal strength")*1e16,
						particle_parameters.at("Bond shear strength")*1e16,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"), particle_parameters.at("Max bond distance"));
	
	eflow.install_wall_bonds(particle_parameters.at("Wall bond normal stiffness")*1e10, 
						particle_parameters.at("Wall bond shear stiffness")*1e10,
						particle_parameters.at("Wall bond normal strength")*1e10,
						particle_parameters.at("Wall bond shear strength")*1e10,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Domain length y")*0.45, particle_parameters.at("Max wall bond distance"));

	// Create a vector to hold particle information
	std::vector<std::vector<double>> particle_info;
	// For a given particle: x pos, y pos, time, repetition
	std::vector<double> temp_info(4);
	// For particle properties
	const std::vector<std::unique_ptr<Particle>>& particles = eflow.get_particles();

	// Bond info
	std::list<Wall_bond>& bonds = eflow.get_wall_bonds();

	for (size_t ir = 0; ir < 200; ++ir)
	{
		// Lift force
		eflow.add_lift_force(0.1, 1e3, 1e-3);
		
		// Align heading with neighboring cells
		eflow.compute_average_heading(particle_parameters.at("Interaction radius"),
				needs_external_flow);
	
		// Update particle positions
		eflow.propagate_particles(needs_external_flow);
	
		// Compute bond forces, check for failure, add contributions
		eflow.process_bonds();		

		// Compute bond forces, check for failure, add contributions
		eflow.process_wall_bonds();		

		for (auto& bond_i : bonds)
		{
			// Get the ID of the bonded particle
			size_t pIDs = bond_i.get_particle_ID();
	
			// Compute particle-wall distance
			Particle temp_particle_i = *(particles.at(pIDs));
	
			std::vector<double> temp_coords = temp_particle_i.get_position();
			double xi = temp_coords.at(0); 			
			double yi = temp_coords.at(1);
			double xj = bond_i.get_x_coord(); 			
			double yj = bond_i.get_y_coord();
		
			// Distance (x - periodic, y - walls)
			double dx = std::min<double>(std::fabs(xi - xj), 
					particle_parameters.at("Domain length x") - std::fabs(xi - xj));
			double dy = yi - yj;
		
			double dsq = dx*dx + dy*dy;
	
			// Break the bond if distance is too large
			if (std::sqrt(dsq) > particle_parameters.at("Max wall bond distance"))
			{
				std::cerr << "Long bond not removed " << std::sqrt(dsq)<< std::endl;
				return false;			
			}
		}

		eflow.update_positions_and_headings();

		eflow.reinstall_bonds(particle_parameters.at("Bond normal stiffness")*1e14, 
						particle_parameters.at("Bond shear stiffness")*1e14,
						particle_parameters.at("Bond normal strength")*1e16,
						particle_parameters.at("Bond shear strength")*1e16,
						particle_parameters.at("Particle radius"),
						particle_parameters.at("Bond max gap"));

	}
	return true;	
}

*/
