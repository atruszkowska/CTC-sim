#include "../../include/geometry/geometry_objects/circle.h"
#include "../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the Box class
 *
 *****************************************************/

// Tests
bool setup_test();

// Supporting functions
// Makes a circle and compares the resulting object with expectations
bool check_a_circle(const double radius, 
					const std::vector<double>& center, 
					const double expected_radius, 
					const std::vector<double>& expected_center);

int main()
{
    test_pass(setup_test(), "Make a circle");
}

// Construct a circle and compare with expected
// There are three circle types considered
bool setup_test()
{
	// Input data
	const std::vector<double> all_dimensions{10.2, 2.1e-7, 1.6};		
	const std::vector<std::vector<double>> all_centers{{0.1, 0.1}, 
										{6.4, 10.5}, {-1.0, 0.0}};
	const std::vector<double> expected_dimensions{10.2, 2.1e-7, 1.6};		
	const std::vector<std::vector<double>> expected_centers{{0.1, 0.1}, 
										{6.4, 10.5}, {-1.0, 0.0}};

	// Run a test on each set
	for (size_t i = 0; i < all_dimensions.size(); ++i) {
		if (!check_a_circle(all_dimensions.at(i), all_centers.at(i),
			expected_dimensions.at(i), expected_centers.at(i))) {
			std::cerr << "Circle making test failed for dataset " << i << std::endl;
			return false;
		}
	}

	return true;
}

bool check_a_circle(const double radius, 
					const std::vector<double>& center, 
					const double expected_radius, 
					const std::vector<double>& expected_center)
{
	double tol = 1e-5;
	Circle test_circle;

	// Dimensions check
	test_circle.set_radius(radius);
	if (!(float_equality(test_circle.R(), expected_radius, tol))) {
		std::cerr << "Circle radius does not match expected" << std::endl;
		return false;
	}

	// Setup the center coordinates
	test_circle.set_center(center.at(0), center.at(1));

	// Center check
	const std::vector<double> cur_center = test_circle.get_center_coordinates();
	if (!std::equal(cur_center.cbegin(), cur_center.cend(), expected_center.cbegin(), 
			[tol](double comp, double exp){ return float_equality(comp, exp, tol); })) {
		std::cerr << "Circle center not equal to expected"<< std::endl;
		return false;
	}
	
	return true;
}
