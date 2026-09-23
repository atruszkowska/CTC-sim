#include "../../include/geometry/geometry_objects/box.h"
#include "../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the Box class
 *
 *****************************************************/

// Tests
bool setup_test();

// Supporting functions
// Makes the box and compares the resulting object with expectations
bool check_a_box(const std::vector<double>&, const std::vector<double>&, 
					const std::vector<double>&, const std::vector<double>&);

int main()
{
    test_pass(setup_test(), "Make a box");
}

// Construct a box and compare with expected
// There are three box types considered
bool setup_test()
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
		if (!check_a_box(all_dimensions.at(i), all_centers.at(i),
			expected_dimensions.at(i), expected_bounds.at(i))) {
			std::cerr << "Box making test failed for dataset " << i << std::endl;
			return false;
		}
	}

	return true;
}

bool check_a_box(const std::vector<double>& dimensions, 
					const std::vector<double>& center, 
					const std::vector<double>& expected_dimensions, 
					const std::vector<double>& expected_bounds)
{
	double tol = 1e-5;
	Box test_box;

	// Dimensions check
	test_box.set_dimensions(dimensions.at(0), dimensions.at(1));
	if (!(float_equality(test_box.Lx(), expected_dimensions.at(0), tol))) {
		std::cerr << "x dimension does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(test_box.Ly(), expected_dimensions.at(1), tol))) {
		std::cerr << "y dimension does not match expected"<< std::endl;
		return false;
	}

	// Setup the center coordinates
	test_box.set_center(center.at(0), center.at(1));

	// Bounds check
	test_box.compute_bounds();
	const std::vector<double> bounds = test_box.get_limits();
	if (!std::equal(bounds.cbegin(), bounds.cend(), expected_bounds.cbegin(), 
			[tol](double comp, double exp){ return float_equality(comp, exp, tol); })) {
		std::cerr << "Bounds not equal to expected"<< std::endl;
		return false;
	}
	
	return true;
}
