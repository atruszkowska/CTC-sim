#include "../../include/boundaries/periodic_boundary.h"
#include "../../include/boundaries/wall_boundary.h"
#include "../../include/boundaries/lj_wall.h"
#include "../../include/io_operations/easy_io.h"
#include "../../include/test_utils.h"
#include "../../include/utils.h"

/***************************************************** 
 *
 * Test suite for the Boundary class hierarchy 
 *
 *****************************************************/

// Tests
bool periodic_boundaries_test();
bool wall_boundaries_test();
bool lj_force_test();
bool lj_boundary_wall_1();
bool lj_boundary_wall_2();
bool lj_boundary_walls_1_and_2();

// Supporting functions
// Make an LJ potential that is fixed besides the wall distance L
// Apply it to particle position x 
void make_and_apply_fixed_LJ(double L, double& x); 

int main()
{
	test_pass(periodic_boundaries_test(), "Periodic boundaries");
	test_pass(wall_boundaries_test(), "Wall boundaries");
	test_pass(lj_force_test(), "LJ force values");
	test_pass(lj_boundary_wall_1(), "LJ Boundary: cell near wall 1, far from wall 2");
	test_pass(lj_boundary_wall_2(), "LJ Boundary: cell near wall 2, far from wall 1");
	test_pass(lj_boundary_walls_1_and_2(), "LJ Boundary: walls 1 and 2 close, cell everywhere");
}

// Verify if PBCs are working correctly
bool periodic_boundaries_test()
{
	double tol = 1e-5;
	std::vector<std::vector<double>> limits = {{0, 1}, {-2, -1}, {1, 3}};
	std::vector<double> lengths = {1, 1, 2};
	std::vector<std::vector<double>> positions = {{-50.7, -0.5, 0.2, 1.3}, 
			{-2.5, -1.8, -0.3, 1.2}, {0.7, 2.5, 3.1, 5.4}};
	std::vector<std::vector<double>> expected = {{0.3, 0.5, 0.2, 0.3}, 
			{-1.5, -1.8, -1.3, -1.8}, {2.7, 2.5, 1.1, 1.4}};

	for (size_t i = 0; i < lengths.size(); ++i) {
		Periodic_boundary pbc(lengths.at(i), limits.at(i));
		for (size_t j = 0; j < 3; ++j) {
			pbc.apply_boundary_condition(positions.at(i).at(j));
			if (!(float_equality(positions.at(i).at(j), expected.at(i).at(j), tol))) {
				std::cerr << "position " << positions.at(i).at(j)  
						  << " does not match expected " << expected.at(i).at(j) 
						  << " after applying the boundary conditions" << std::endl;
				return false;
			}					
		}	
	}

	return true;
}

// Verify if wall boundaries are working correctly
bool wall_boundaries_test()
{
	double tol = 1e-5;
	std::vector<std::vector<double>> limits = {{0, 1}, {-2, -1}, {1, 3}};
	std::vector<double> lengths = {1, 1, 2};
	std::vector<double> R_min = {0.1, 0.01, 1.0};
	std::vector<std::vector<double>> positions = {{-50.7, -0.5, 0.2, 1.3}, 
			{-2.5, -1.8, -0.3, 1.2}, {0.7, 2.5, 3.1, 5.4}};
	std::vector<std::vector<double>> expected = {{0.1, 0.1, 0.2, 0.9}, 
			{-1.99, -1.8, -1.01, -1.01}, {2.0, 2.0, 2.0, 2.0}};

	for (size_t i = 0; i < lengths.size(); ++i) {
		Wall_boundary wbc(lengths.at(i), limits.at(i), R_min.at(i));
		for (size_t j = 0; j < 3; ++j) {
			wbc.apply_boundary_condition(positions.at(i).at(j));
			if (!(float_equality(positions.at(i).at(j), expected.at(i).at(j), tol))) {
				std::cerr << "position " << positions.at(i).at(j)  
						  << " does not match expected " << expected.at(i).at(j) 
						  << " after applying the boundary conditions" << std::endl;
				return false;
			}					
		}	
	}

	return true;
}

// Verify that the LJ force is correctly computed
bool lj_force_test()
{
	const double L = 200e-6, epsilon = 0.8, alpha = 2, sigma = 1e-5, dt = 0.1, m_cell = 0.1; 
	const double tol = 1e-8;
	const std::vector<double> lim{0, L};
	const int m = 9, n = 3; 

	// C++ LJ
	Lj_wall ljw(L, lim, sigma, epsilon, alpha, m, n, dt, m_cell);
	const std::vector<double> lj_force = ljw.compute_lj_force(1e-8, L, 100);

	// Actual solution
	Easy_IO io("lj_force_data.txt");
	std::vector<std::vector<double>> temp_solution = io.read_vector<double>();
	std::vector<double> true_solution = temp_solution.at(0);

	// Comparison	
	bool are_correct = is_equal_floats<double>({lj_force}, {true_solution}, tol);
	if (!are_correct) 
	{
		std::cerr << "Computed LJ force does not match expected" << std::endl;
	    for (size_t i = 0; i < lj_force.size(); ++i) {
	        std::cout <<  lj_force.at(i) << " " << true_solution.at(i) << "\n";
	    }
		return are_correct;
	}
	return are_correct;
}

// Suite 1: walls far apart, impact of Wall 1 
bool lj_boundary_wall_1()
{
	// Initial postions (changed by BCs)
	std::vector<double> x0 = {1.01e-7, 2.99e-7, 2.0099e-5};
	// Final positions
	const std::vector<double> xf = {395811, -2.2404e2, 2.0099e-5};
	// Wall distance
	const double L = 100;
	// Comparison tolerance
	const double tol = 1e-5;

	for (size_t i = 0; i < x0.size(); ++i)
	{
		make_and_apply_fixed_LJ(L, x0.at(i));
		if (!float_equality(xf.at(i), x0.at(i), tol)) 
		{
			std::cerr << "Expected " << xf.at(i) << " got " << x0.at(i) << "\n";
	        	return false;
	    }
	}
	return true;
}

// Suite 2: walls far apart, impact of Wall 2 
bool lj_boundary_wall_2()
{
	// Wall distance
	const double L = 100;
	// Initial postions (changed by BCs)
	std::vector<double> x0 = {L - 1.01e-7, L - 2.99e-7, L - 2.0099e-5};
	// Final positions
	const std::vector<double> xf = {-395711, 3.2404e2, 100};
	// Comparison tolerance
	const double tol = 1e-5;

	for (size_t i = 0; i < x0.size(); ++i)
	{
		make_and_apply_fixed_LJ(L, x0.at(i));
		if (!float_equality(xf.at(i), x0.at(i), tol)) 
		{
			std::cerr << "Expected " << xf.at(i) << " got " << x0.at(i) << "\n";
	        	return false;
	    }
	}
	return true;
}

// Suite 3: walls close, impact of both walls
bool lj_boundary_walls_1_and_2()
{
	// Wall distance
	const double L = 1e-6;
	// Initial postions (changed by BCs)
	std::vector<double> x0 = {1.01e-7, 2.99e-7, 5.0099e-7, L - 1.01e-7, L - 2.99e-7};
	// Final positions
	const std::vector<double> xf = {3.9581e+05, -2.23463e+02, 1.7028e-01, -3.9581e+05, 2.23463e+02};
	// Comparison tolerance
	const double tol = 1e-5;

	for (size_t i = 0; i < x0.size(); ++i)
	{
		make_and_apply_fixed_LJ(L, x0.at(i));
		if (!float_equality(xf.at(i), x0.at(i), tol)) 
		{
			std::cerr << "Expected " << xf.at(i) << " got " << x0.at(i) << "\n";
	        	return false;
	    }
	}
	return true;
}

// Make an LJ potential that is fixed besides the wall distance L
// Apply it to particle position x
void make_and_apply_fixed_LJ(double L, double& x)
{
	const double epsilon = 0.8, alpha = 0.1, sigma = 1e-7, dt = 0.1, m_cell = 0.1; 
	const std::vector<double> lim{0, L};
	const int m = 12, n = 6; 

	Lj_wall ljw(L, lim, sigma, epsilon, alpha, m, n, dt, m_cell);
	ljw.apply_boundary_condition(x);
}

