#include "../../include/geometry/geometry.h"
#include "../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the Geometry class
 *
 *****************************************************/

// Tests
bool setup_test();
bool well_setup_test();

// Supporting functions
// Makes the box and compares the resulting object with expectation
bool check_a_box(const std::vector<double>&, const std::vector<double>&, 
					const std::vector<double>&, const std::vector<double>&);
// Supporting functions
// Makes a circle and compares the resulting object with expectations
bool check_a_circle(const double radius, 
					const std::vector<double>& center, 
					const double expected_radius, 
					const std::vector<double>& expected_center);

int main()
{
    test_pass(setup_test(), "Make a rectangular domain");
    test_pass(well_setup_test(), "Make a circular domain");
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

// Construct a well (circle) and compare with expected
// There are three circle types considered
bool well_setup_test()
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

bool check_a_box(const std::vector<double>& dimensions,
					const std::vector<double>& center,	
					const std::vector<double>& expected_dimensions, 
					const std::vector<double>& expected_bounds)
{
	double tol = 1e-5;
	Geometry test_geometry;

	// Dimensions check
	test_geometry.make_box(dimensions.at(0), dimensions.at(1), center.at(0), center.at(1));
	if (!(float_equality(test_geometry.get_x_length(), expected_dimensions.at(0), tol))) {
		std::cerr << "x dimension does not match expected" << std::endl;
		return false;
	}
	if (!(float_equality(test_geometry.get_y_length(), expected_dimensions.at(1), tol))) {
		std::cerr << "y dimension does not match expected"<< std::endl;
		return false;
	}


	// Bounds check
	const std::vector<double> bounds = test_geometry.get_box_limits();
	if (!std::equal(bounds.cbegin(), bounds.cend(), expected_bounds.cbegin(), 
			[tol](double comp, double exp){ return float_equality(comp, exp, tol); })) {
		std::cerr << "Bounds not equal to expected"<< std::endl;
		return false;
	}
	
	return true;
}

bool check_a_circle(const double radius, 
					const std::vector<double>& center, 
					const double expected_radius, 
					const std::vector<double>& expected_center)
{
	double tol = 1e-5;
	Geometry geom;

	geom.make_circle(radius, center.at(0), center.at(1));

	// Dimensions check
	if (!(float_equality(geom.get_radius(), expected_radius, tol))) {
		std::cerr << "Circle radius does not match expected" << std::endl;
		return false;
	}

	// Center check
	const std::vector<double> cur_center = geom.get_center();
	if (!std::equal(cur_center.cbegin(), cur_center.cend(), expected_center.cbegin(), 
			[tol](double comp, double exp){ return float_equality(comp, exp, tol); })) {
		std::cerr << "Circle center not equal to expected"<< std::endl;
		return false;
	}
	
	return true;
}
