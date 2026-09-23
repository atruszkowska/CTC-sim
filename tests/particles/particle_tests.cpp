#include "../../include/particles/particle.h"
#include "../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the Particle class
 *
 *****************************************************/

// Tests
bool constructor_tests();
bool check_velocity();
bool check_new_position();
bool check_heading();


// Supporting functions
bool make_and_check(const double x0, const double y0, 
				const double theta0, const double vmag);

int main()
{
    test_pass(constructor_tests(), "Make a particle");
    test_pass(check_velocity(), "Check particle velocity");
    test_pass(check_new_position(), "Check particle position");
    test_pass(check_heading(), "Check particle heading");
}

// Create a particle with and without velocity
bool constructor_tests()
{
	// Stationary
	if (!make_and_check(1.1, 1.2, 0.2, 0.0)) {
		std::cerr << "Test for the stationary particle failed" << std::endl;
		return false;
	}
		
	// Moving
	if (!make_and_check(0.3, 10.7, 0.5, 0.25)) {
		std::cerr << "Test for the moving particle failed" << std::endl;
		return false;
	}

	return true;
}

// Construct a particle and check if its properties are correct
bool make_and_check(const double x0, const double y0, 
				const double theta0, const double vmag)
{
	double tol = 1e-5;

	Particle test_particle(x0, y0, theta0, vmag);
	
	const std::vector<double> pos = test_particle.get_position();
	if (!(float_equality(pos.at(0), x0, tol))) {
		std::cerr << "x coordinate does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(pos.at(1), y0, tol))) {
		std::cerr << "y coordinate does not match expected" << std::endl;
		return false;
	}

	const std::vector<double> vel = test_particle.get_velocity();
	if (!(float_equality(vel.at(0), vmag*std::cos(theta0), tol))) {
		std::cerr << "x velocity component does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(vel.at(1), vmag*std::sin(theta0), tol))) {
		std::cerr << "y velocity component does not match expected" << std::endl;
		return false;
	}
	
	const double pos_x = test_particle.get_x_position();
	const double pos_y = test_particle.get_y_position();
	if (!(float_equality(pos_x, x0, tol))) {
		std::cerr << "x coordinate does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(pos_y, y0, tol))) {
		std::cerr << "y coordinate does not match expected" << std::endl;
		return false;
	}

	const double vel_mag = test_particle.get_velocity_magnitude();
	if (!(float_equality(vel_mag, vmag, tol))) {
		std::cerr << "velocity magnitude does not match expected" << std::endl;
		return false;
	}


	const double heading = test_particle.get_heading();
	if (!(float_equality(heading, theta0, tol))) {
		std::cerr << "Heading does not match expected" << std::endl;
		return false;
	}

	return true;
}

// Dynamic properties of particle
bool check_velocity() {
	
	double tol = 1e-5;
	double test_vel_mag = 0.25;
	double test_heading = 0.5;
	double test_pos_x = 0.3;
	double test_pos_y = 10.7;

	Particle t_particle(test_pos_x, test_pos_y, test_heading, test_vel_mag);
	t_particle.compute_velocity();

	// Access velocity components
	const std::vector<double> vel = t_particle.get_velocity();

	// Define expected values
	const double expected_vel_x = test_vel_mag * std::cos(test_heading);
	const double expected_vel_y = test_vel_mag * std::sin(test_heading);

	// Check computed velocity components
	if (!(float_equality(vel.at(0), expected_vel_x, tol))) {
		std::cerr << "x component of the velocity does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(vel.at(1), expected_vel_y, tol))) {
		std::cerr << "y component of the velocity does not match expected" << std::endl;
		return false;
	}
	return true;
}

bool check_new_position() {
	
	double tol = 1e-5;
	Particle t_particle(0.3, 10.7, 0.5, 0.25);

	double test_dt = 2.0;
	double test_vel_mag = 0.25;
	double test_heading = 0.5;
	double test_pos_x = 0.3;
	double test_pos_y = 10.7;

	
	// Access particle position
	t_particle.compute_new_position(test_dt);
	const std::vector<double> pos = t_particle.get_position();

	// Define expected values
	const double expected_x = test_pos_x + test_vel_mag * std::cos(test_heading) * test_dt;
	const double expected_y = test_pos_y + test_vel_mag * std::sin(test_heading) * test_dt;

	// Check computed particle position
	if (!(float_equality(pos.at(0), expected_x, tol))) {
		std::cerr << "x coordinate does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(pos.at(1), expected_y, tol))) {
		std::cerr << "y coordinate does not match expected" << std::endl;
		return false;
	}
	return true;
}

bool check_heading() {
	
	double tol = 1e-5;
	double theta_1 = 0.9;
	const double noise = 1.0;
	const std::vector<double> new_pos = {1.5, 3.2};
	const double expected_heading = theta_1 + noise;
	
	double test_vel_mag = 0.25;
	double test_heading = 0.5;
	double test_pos_x = 0.3;
	double test_pos_y = 10.7;
	Particle t_particle(test_pos_x, test_pos_y, test_heading, test_vel_mag);

	t_particle.set_new_heading(std::move(theta_1));
	t_particle.update_heading();
	t_particle.add_noise(noise);
	t_particle.set_new_position(std::move(new_pos));

	const double heading = t_particle.get_heading();
	if (!(float_equality(heading, expected_heading, tol))) {
		std::cerr << "Incorrect after setting and updating heading" << std::endl;
		return false;
	}

	const std::vector<double> pos = t_particle.get_position();
	if (!(float_equality(pos.at(0), new_pos.at(0), tol))) {
		std::cerr << "Incorrect value after setting new x coordinate" << std::endl;
		return false;
	}
	if (!(float_equality(pos.at(1), new_pos.at(1), tol))) {
		std::cerr << "Incorrect value after setting new y coordinate" << std::endl;
		return false;
	}

	return true;
}

