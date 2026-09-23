#include "../../include/interactions/basic.h"
#include "../../include/interactions/bond.h"
#include "../../include/interactions/wall_bond.h"
#include "../../include/geometry/geometry.h"
#include "../../include/particles/particle.h"
#include "../../include/rng.h"
#include "../../include/test_utils.h"

// Tests
bool basic_overlap_test();
bool bond_test();
bool bond_failure_test();
bool bond_contribution_test();
bool wall_bond_test();
bool wall_bond_failure_test();
bool wall_bond_contribution_test();

// Supporting functions
double compute_particle_distance(const Particle&, const Particle&);
double compute_particle_distance(const Particle& pt1, const std::vector<double> pos_2);

int main()
{
	test_pass(basic_overlap_test(), "Basic removal of overlaps");
	test_pass(bond_test(), "Bond formation and application");
	test_pass(bond_failure_test(), "Bond failure");
	test_pass(bond_contribution_test(), "Bond contribution to particle motion");
	test_pass(wall_bond_test(), "Wall bond formation and application");
	test_pass(wall_bond_failure_test(), "Wall bond failure");
	test_pass(wall_bond_contribution_test(), "Wall bond contribution to particle motion");
}

bool basic_overlap_test()
{
	// For floating point comparisons
	const double tol = 1e-5;

	// Random number generator
	RNG rng;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// Interactions
	// Minimum distance between particles
	double particle_d_min = 0.75, distance = -1.0;
	Basic interactions;	

	// Two overlaping particles
	const double x1 = 1.0, y1 = 2.0, x2 = 1.25, y2 = 2.5; 
	const double theta1 = 0.1, theta2 = 0.0, vmag = 1.5;
	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, theta1, vmag)));
	particles.push_back(std::unique_ptr<Particle>(new Particle(x2, y2, theta2, vmag)));

	// First check if really overlapping, then remove the overlap and check again
	distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	if ( distance >= particle_d_min ) 
	{
		std::cerr << "The distance between particles is not smaller than minimum allowable " 
				  << distance << std::endl;
		return false;
	}
	interactions.remove_overlaps(particles.at(0), particles.at(1), rng, distance*distance, particle_d_min);

	distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	if ( distance < particle_d_min ) 
	{
		std::cerr << "The overlap between particles was not removed " 
				  << distance << " (Required minimum: " << particle_d_min << ")"<< std::endl;
		return false;
	}

	// No change now
	interactions.remove_overlaps(particles.at(0), particles.at(1), rng, distance*distance, particle_d_min);
	double new_distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	if (!float_equality(new_distance, distance, tol)) 
	{
		std::cerr << "The distance between particles 1 and 2 shouldn't change" << std::endl;
		return false;
	}

	// Move 2 so it's at the exact same spot as 1
	particles.at(1)->set_new_position({particles.at(0)->get_x_position(),
					particles.at(0)->get_y_position()});
	distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	if ( distance >= particle_d_min ) 
	{
		std::cerr << "The distance between particles is not smaller than minimum allowable " 
				  << distance << std::endl;
		return false;
	}		
	interactions.remove_overlaps(particles.at(0), particles.at(1), rng, distance*distance, particle_d_min);
	distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	if ( distance < particle_d_min ) 
	{
		std::cerr << "The overlap between particles was not removed " 
				  << distance << " (Required minimum: " << particle_d_min << ")"<< std::endl;
		return false;
	}
	if ( std::isnan(distance) ) 
	{
		std::cerr << "One particle position is now NaN " 
				  << distance << std::endl;
		return false;
	}
	return true;
}

bool bond_test()
{
	// For floating point comparisons
	const double tol = 1e-5;

	// Random number generator
	RNG rng;

	// Time step
	const double delta_t = 10.0;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// Two overlaping particles with angular and translational velocities
	const double x1 = 0.1, y1 = 1.5, x2 = 0.5, y2 = 1.8; 
	const double omega_1 = 1.5, omega_2 = 3.0;
	const std::vector<double> vel_1{0.25, 0.1}, vel_2{0.5, 0.3};
	const double R_particle = 2.0;

	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, 0, 0)));
	particles.push_back(std::unique_ptr<Particle>(new Particle(x2, y2, 0, 0)));
	
	// Set angular and translational velocities
	particles.at(0)->set_angular_velocity(omega_1);
	particles.at(1)->set_angular_velocity(omega_2);
	particles.at(0)->set_translational_velocity(vel_1);
	particles.at(1)->set_translational_velocity(vel_2);

	// Bond properties
	const double k_n = 100;
	const double k_s = 10;
	const double sigma_max = 2500;
	const double tau_max = 150;

	// Install bond
	Bond particle_bond(k_n, k_s, sigma_max, tau_max, R_particle, 0, 1);

	// Compute force and moment increments
	const double distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	particle_bond.add_displacement_increments(particles.at(0), particles.at(1), rng, distance, delta_t);

	// Expected values
	double F_n_exp = 1280;
	std::vector<double> F_s_exp{-267.6, 356.8};
	double M_b_exp = -8000.0;

	// Compare with expected 
	double F_n = particle_bond.get_normal_force();
	std::vector<double> F_s = particle_bond.get_shear_force();
	double M_b = particle_bond.get_bending_moment();

	if (!float_equality(F_n, F_n_exp, tol)) 
	{
		std::cerr << "Normal force " << F_n << " not equal to expected " << F_n_exp << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(0), F_s_exp.at(0), tol)) 
	{
		std::cerr << "x component of shear force " << F_s.at(0) << " not equal to expected " << F_s_exp.at(0) << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(1), F_s_exp.at(1), tol)) 
	{
		std::cerr << "y component of shear force " << F_s.at(1) << " not equal to expected " << F_s_exp.at(1) << std::endl;
		return false;
	}
	if (!float_equality(M_b, M_b_exp, tol)) 
	{
		std::cerr << "Bending moment " << M_b << " not equal to expected " << M_b_exp << std::endl;
		return false;
	}

	// Add one more increment, compare again
	particle_bond.add_displacement_increments(particles.at(0), particles.at(1), rng, distance, delta_t);

	// Expected values
	F_n_exp = 2*1280;
	F_s_exp = {-267.6*2, 356.8*2};
	M_b_exp = -8000.0*2;

	// Compare with expected 
	F_n = particle_bond.get_normal_force();
	F_s = particle_bond.get_shear_force();
	M_b = particle_bond.get_bending_moment();

	if (!float_equality(F_n, F_n_exp, tol)) 
	{
		std::cerr << "Second increment: normal force " << F_n << " not equal to expected " << F_n_exp << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(0), F_s_exp.at(0), tol)) 
	{
		std::cerr << "Second increment: x component of shear force " << F_s.at(0) << " not equal to expected " << F_s_exp.at(0) << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(1), F_s_exp.at(1), tol)) 
	{
		std::cerr << "Second increment: y component of shear force " << F_s.at(1) << " not equal to expected " << F_s_exp.at(1) << std::endl;
		return false;
	}
	if (!float_equality(M_b, M_b_exp, tol)) 
	{
		std::cerr << "Second increment: bending moment " << M_b << " not equal to expected " << M_b_exp << std::endl;
		return false;
	}

	return true;
}

bool bond_failure_test()
{
	// Random number generator
	RNG rng;

	// Time step
	const double delta_t = 10.0;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// Two overlaping particles with angular and translational velocities
	const double x1 = 0.1, y1 = 1.5, x2 = 0.5, y2 = 1.8; 
	const double omega_1 = 1.5, omega_2 = 3.0;
	const std::vector<double> vel_1{0.25, 0.1}, vel_2{0.5, 0.3};
	const double R_particle = 2.0;

	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, 0, 0)));
	particles.push_back(std::unique_ptr<Particle>(new Particle(x2, y2, 0, 0)));
	
	// Set angular and translational velocities
	particles.at(0)->set_angular_velocity(omega_1);
	particles.at(1)->set_angular_velocity(omega_2);
	particles.at(0)->set_translational_velocity(vel_1);
	particles.at(1)->set_translational_velocity(vel_2);

	// Bond properties
	const double k_n = 100;
	const double k_s = 10;
	double sigma_max = 2700;
	double tau_max = 150;

	// Install bond
	Bond particle_bond(k_n, k_s, sigma_max, tau_max, R_particle, 0, 1);

	// Compute force and moment increments
	const double distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	particle_bond.add_displacement_increments(particles.at(0), particles.at(1), rng, distance, delta_t);

	// Failure variables
	bool failed = false, failed_in_tension = false, failed_in_shear = false;

	// This bond shouldn't fail
	failed = particle_bond.bond_failed(failed_in_shear, failed_in_tension);
	if (failed || failed_in_shear || failed_in_tension)
	{
		std::cerr << "The bond should not have failed" << std::endl;
		return false;
	}

	// This bond should fail in tension
	// Install bond
	sigma_max = 2600;
	Bond particle_bond_tensile_failure(k_n, k_s, sigma_max, tau_max, R_particle, 0, 1);
	// Compute force and moment increments
	particle_bond_tensile_failure.add_displacement_increments(particles.at(0), particles.at(1), rng, distance, delta_t);
	failed = particle_bond_tensile_failure.bond_failed(failed_in_shear, failed_in_tension);
	if (!failed || failed_in_shear || !failed_in_tension)
	{
		std::cerr << "This bond should failed in tension" << std::endl;
		return false;
	}
	failed = false;
	failed_in_shear = false;
	failed_in_tension = false; 

	// This bond should fail in tension
	// Install bond
	sigma_max = 3000;
	tau_max = 100;
	Bond particle_bond_shear_failure(k_n, k_s, sigma_max, tau_max, R_particle, 0, 1);
	// Compute force and moment increments
	particle_bond_shear_failure.add_displacement_increments(particles.at(0), particles.at(1), rng, distance, delta_t);
	failed = particle_bond_shear_failure.bond_failed(failed_in_shear, failed_in_tension);
	if (!failed || !failed_in_shear || failed_in_tension)
	{
		std::cerr << "This bond should failed in shear" << std::endl;
		return false;
	}
	failed = false;
	failed_in_shear = false;
	failed_in_tension = false;

	return true;
}

// Tests if bond forces and moment are correctly accounted for in particle motion
bool bond_contribution_test()
{
	// For floating point comparisons
	const double tol = 1e-5;

	// Random number generator
	RNG rng;

	// Time step
	const double delta_t = 10.0;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// Two overlaping particles with angular and translational velocities
	const double x1 = 0.1, y1 = 1.5, x2 = 0.5, y2 = 1.8, mass = 700; 
	const double omega_1 = 1.5, omega_2 = 3.0;
	const double theta_1 = 3.7, theta_2 = 1.2;
	const std::vector<double> vel_1{0.25, 0.1}, vel_2{0.5, 0.3};
	const double R_particle = 2.0;

	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, 0, 0)));
	particles.push_back(std::unique_ptr<Particle>(new Particle(x2, y2, 0, 0)));
	
	// Set angular and translational velocities
	particles.at(0)->set_angular_velocity(omega_1);
	particles.at(1)->set_angular_velocity(omega_2);
	particles.at(0)->set_translational_velocity(vel_1);
	particles.at(1)->set_translational_velocity(vel_2);
	particles.at(0)->set_heading(theta_1);
	particles.at(1)->set_heading(theta_2);
	particles.at(0)->set_new_heading(theta_1);
	particles.at(1)->set_new_heading(theta_2);
	particles.at(0)->set_particle_mass(mass);
	particles.at(1)->set_particle_mass(mass);

	// Bond properties
	const double k_n = 100;
	const double k_s = 10;
	const double sigma_max = 25000;
	const double tau_max = 1500;

	// Install bond
	Bond particle_bond(k_n, k_s, sigma_max, tau_max, R_particle, 0, 1);

	// Compute force and moment increments
	const double distance = compute_particle_distance(*(particles.at(0)),*(particles.at(1)));
	particle_bond.add_displacement_increments(particles.at(0), particles.at(1), rng, distance, delta_t);

	// Include bond influence in particle position and heading
	particle_bond.correct_particle_positions_and_headings(particles.at(0), particles.at(1), delta_t);

	// Update actual positions and headings
	particles.at(0)->update_heading();
	particles.at(1)->update_heading();
	particles.at(0)->compute_new_position();
	particles.at(1)->compute_new_position();

	// Expected values
	const std::vector<double> pos_j_exp{-107.5571, -158.8857}, pos_i_exp{108.1571, 162.1857};
	const double theta_j_exp = 572.6286, theta_i_exp = -567.7286;

	// Actual values
	const std::vector<double> pos_i = particles.at(0)->get_position();
	const std::vector<double> pos_j = particles.at(1)->get_position();
	const double theta_i = particles.at(0)->get_heading();
	const double theta_j = particles.at(1)->get_heading();

	// Compare with expected 
	if (!is_equal_floats<double>({pos_i}, {pos_i_exp}, tol)) 
	{
		std::cerr << "Position of particle i " << pos_i.at(0) << " " << 
			pos_i.at(1) << " not equal to expected " << pos_i_exp.at(0) << " " 
			<< pos_i_exp.at(1) << std::endl;
		return false;
	}
	if (!is_equal_floats<double>({pos_j}, {pos_j_exp}, tol)) 
	{
		std::cerr << "Position of particle j " << pos_j.at(0) << " " << 
			pos_j.at(1) << " not equal to expected " << pos_j_exp.at(0) << " " 
			<< pos_j_exp.at(1) << std::endl;
		return false;
	}
	if (!float_equality(theta_i, theta_i_exp, tol)) 
	{
		std::cerr << "Heading of particle i " << theta_i << " not equal to expected " << theta_i_exp << std::endl;
		return false;
	}
	if (!float_equality(theta_j, theta_j_exp, tol)) 
	{
		std::cerr << "Heading of particle j " << theta_j << " not equal to expected " << theta_j_exp << std::endl;
		return false;
	}

	return true;
}

// Wall bonds
bool wall_bond_test()
{
	// For floating point comparisons
	const double tol = 1e-5;

	// Random number generator
	RNG rng;

	// Time step
	const double delta_t = 10.0;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// One particle with angular and translational velocities
	const double x1 = 0.1, y1 = 1.5, xw = 0.5, yw = 1.8; 
	const double omega_1 = 1.5;
	const std::vector<double> vel_1{0.25, 0.1};
	const double R_particle = 2.0;

	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, 0, 0)));
	
	// Set angular and translational velocities
	particles.at(0)->set_angular_velocity(omega_1);
	particles.at(0)->set_translational_velocity(vel_1);

	// Bond properties
	const double k_n = 100;
	const double k_s = 10;
	const double sigma_max = 2500;
	const double tau_max = 150;

	// Install bond
	Wall_bond wall_bond(k_n, k_s, sigma_max, tau_max, R_particle, 0, {xw, yw});

	// Compute force and moment increments
	const double distance = compute_particle_distance(*(particles.at(0)), {xw, yw});
	wall_bond.add_displacement_increments(particles.at(0), rng, distance, delta_t);

	// Expected values
	double F_n_exp = -1040;
	std::vector<double> F_s_exp{-73.2, 97.6};
	double M_b_exp = 8000.0;

	// Compare with expected 
	double F_n = wall_bond.get_normal_force();
	std::vector<double> F_s = wall_bond.get_shear_force();
	double M_b = wall_bond.get_bending_moment();

	if (!float_equality(F_n, F_n_exp, tol)) 
	{
		std::cerr << "Normal force " << F_n << " not equal to expected " << F_n_exp << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(0), F_s_exp.at(0), tol)) 
	{
		std::cerr << "x component of shear force " << F_s.at(0) << " not equal to expected " << F_s_exp.at(0) << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(1), F_s_exp.at(1), tol)) 
	{
		std::cerr << "y component of shear force " << F_s.at(1) << " not equal to expected " << F_s_exp.at(1) << std::endl;
		return false;
	}
	if (!float_equality(M_b, M_b_exp, tol)) 
	{
		std::cerr << "Bending moment " << M_b << " not equal to expected " << M_b_exp << std::endl;
		return false;
	}

	// Add one more increment, compare again
	wall_bond.add_displacement_increments(particles.at(0), rng, distance, delta_t);

	// Expected values
	F_n_exp *= 2;
	F_s_exp = {-73.2*2, 97.6*2};
	M_b_exp = 8000.0*2;

	// Compare with expected 
	F_n = wall_bond.get_normal_force();
	F_s = wall_bond.get_shear_force();
	M_b = wall_bond.get_bending_moment();

	if (!float_equality(F_n, F_n_exp, tol)) 
	{
		std::cerr << "Second increment: normal force " << F_n << " not equal to expected " << F_n_exp << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(0), F_s_exp.at(0), tol)) 
	{
		std::cerr << "Second increment: x component of shear force " << F_s.at(0) << " not equal to expected " << F_s_exp.at(0) << std::endl;
		return false;
	}
	if (!float_equality(F_s.at(1), F_s_exp.at(1), tol)) 
	{
		std::cerr << "Second increment: y component of shear force " << F_s.at(1) << " not equal to expected " << F_s_exp.at(1) << std::endl;
		return false;
	}
	if (!float_equality(M_b, M_b_exp, tol)) 
	{
		std::cerr << "Second increment: bending moment " << M_b << " not equal to expected " << M_b_exp << std::endl;
		return false;
	}

	return true;
}

bool wall_bond_failure_test()
{
	// Random number generator
	RNG rng;

	// Time step
	const double delta_t = 10.0;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// One particle with angular and translational velocities
	const double x1 = 0.1, y1 = 1.5, xw = 0.5, yw = 1.8; 
	const double omega_1 = 1.5;
	const std::vector<double> vel_1{0.25, 0.1};
	const double R_particle = 2.0;

	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, 0, 0)));
	
	// Set angular and translational velocities
	particles.at(0)->set_angular_velocity(omega_1);
	particles.at(0)->set_translational_velocity(vel_1);

	// Bond properties
	const double k_n = 100;
	const double k_s = 10;
	double sigma_max = 3261;
	double tau_max = 31;

	// Install bond
	Wall_bond wall_bond(k_n, k_s, sigma_max, tau_max, R_particle, 0, {xw, yw});

	// Compute force and moment increments
	const double distance = compute_particle_distance(*(particles.at(0)), {xw, yw});
	wall_bond.add_displacement_increments(particles.at(0), rng, distance, delta_t);

	// Failure variables
	bool failed = false, failed_in_tension = false, failed_in_shear = false;

	// This bond shouldn't fail
	failed = wall_bond.bond_failed(failed_in_shear, failed_in_tension);
	if (failed || failed_in_shear || failed_in_tension)
	{
		std::cerr << "The bond should not have failed" << std::endl;
		return false;
	}

	// This bond should fail in tension
	// Install bond
	sigma_max = 3260;
	Wall_bond wall_bond_tensile_failure(k_n, k_s, sigma_max, tau_max, R_particle, 0, {xw, yw});
	// Compute force and moment increments
	wall_bond_tensile_failure.add_displacement_increments(particles.at(0), rng, distance, delta_t);
	failed = wall_bond_tensile_failure.bond_failed(failed_in_shear, failed_in_tension);
	if (!failed || failed_in_shear || !failed_in_tension)
	{
		std::cerr << "This bond should failed in tension" << std::endl;
		return false;
	}
	failed = false;
	failed_in_shear = false;
	failed_in_tension = false; 

	// This bond should fail in tension
	// Install bond
	sigma_max = 30000;
	tau_max = 30;
	Wall_bond wall_bond_shear_failure(k_n, k_s, sigma_max, tau_max, R_particle, 0, {xw, yw});
	// Compute force and moment increments
	wall_bond_shear_failure.add_displacement_increments(particles.at(0), rng, distance, delta_t);
	failed = wall_bond_shear_failure.bond_failed(failed_in_shear, failed_in_tension);
	if (!failed || !failed_in_shear || failed_in_tension)
	{
		std::cerr << "This bond should failed in shear" << std::endl;
		return false;
	}
	failed = false;
	failed_in_shear = false;
	failed_in_tension = false;

	return true;
}

// Wall bonds
bool wall_bond_contribution_test()
{
	// For floating point comparisons
	const double tol = 1e-5;

	// Random number generator
	RNG rng;

	// Time step
	const double delta_t = 10.0;

	// Geometry
	const double Lx = 5.0, Ly = 10.0, xc = 0.0, yc = 5.0; 
	Geometry geom;
	geom.make_box(Lx, Ly, xc, yc);

	// Particles 
	std::vector<std::unique_ptr<Particle>> particles;

	// One particle with angular and translational velocities
	const double x1 = 0.1, y1 = 1.5, xw = 0.5, yw = 1.8, mass = 700; 
	const double theta_1 = 3.7, omega_1 = 1.5;
	const std::vector<double> vel_1{0.25, 0.1};
	const double R_particle = 2.0;

	particles.push_back(std::unique_ptr<Particle>(new Particle(x1, y1, 0, 0)));
	
	// Set angular and translational velocities
	particles.at(0)->set_angular_velocity(omega_1);
	particles.at(0)->set_translational_velocity(vel_1);
	particles.at(0)->set_heading(theta_1);
	particles.at(0)->set_new_heading(theta_1);
	particles.at(0)->set_particle_mass(mass);

	// Bond properties
	const double k_n = 100;
	const double k_s = 10;
	const double sigma_max = 2500;
	const double tau_max = 150;

	// Install bond
	Wall_bond wall_bond(k_n, k_s, sigma_max, tau_max, R_particle, 0, {xw, yw});

	// Compute force and moment increments
	const double distance = compute_particle_distance(*(particles.at(0)), {xw, yw});
	wall_bond.add_displacement_increments(particles.at(0), rng, distance, delta_t);

	// Include bond influence in particle position and heading
	wall_bond.correct_particle_positions_and_headings(particles.at(0), delta_t);

	// Update actual positions and headings
	particles.at(0)->update_heading();
	particles.at(0)->compute_new_position();

	// Expected values
	const std::vector<double> pos_i_exp{-129.214, -73.7};
	const double theta_i_exp = 575.1286;

	// Actual values
	const std::vector<double> pos_i = particles.at(0)->get_position();
	const double theta_i = particles.at(0)->get_heading();

	// Compare with expected 
	if (!is_equal_floats<double>({pos_i}, {pos_i_exp}, tol)) 
	{
		std::cerr << "Position of particle i " << pos_i.at(0) << " " << 
			pos_i.at(1) << " not equal to expected " << pos_i_exp.at(0) << " " 
			<< pos_i_exp.at(1) << std::endl;
		return false;
	}
	if (!float_equality(theta_i, theta_i_exp, tol)) 
	{
		std::cerr << "Heading of particle i " << theta_i << " not equal to expected " << theta_i_exp << std::endl;
		return false;
	}

	return true;
}

// Compute distance between two particles
double compute_particle_distance(const Particle& pt1, const Particle& pt2)
{
	const std::vector<double> pos_1 = pt1.get_position();
	const std::vector<double> pos_2 = pt2.get_position();

	return std::sqrt((pos_1.at(0) - pos_2.at(0))*(pos_1.at(0) - pos_2.at(0)) 
		+ (pos_1.at(1) - pos_2.at(1))*(pos_1.at(1) - pos_2.at(1)));
}

// Compute distance between a particle and a wall 
double compute_particle_distance(const Particle& pt1, const std::vector<double> pos_2)
{
	const std::vector<double> pos_1 = pt1.get_position();

	return std::sqrt((pos_1.at(0) - pos_2.at(0))*(pos_1.at(0) - pos_2.at(0)) 
		+ (pos_1.at(1) - pos_2.at(1))*(pos_1.at(1) - pos_2.at(1)));
}
