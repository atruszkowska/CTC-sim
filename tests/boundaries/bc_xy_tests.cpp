#include "../../include/boundaries/round_wall.h"
#include "../../include/io_operations/easy_io.h"
#include "../../include/test_utils.h"
#include "../../include/utils.h"

/***************************************************** 
 *
 * Test suite for the Boundary_xy class hierarchy 
 *
 *****************************************************/

// Tests
bool round_wall_test();

int main()
{
	test_pass(round_wall_test(), "Round enclosure (wall) boundary");
}

// Verify if round wall BCs are working correctly
bool round_wall_test()
{
    double tol = 1e-5;
    double x_c = 1.0;
    double y_c = -2.0;
    double R_min = 0.1;
	double R = 3.0;
	double R_max = R - R_min; 

    Round_wall wall(R, {x_c, y_c}, R_min);

	// Test cases: {x_in, y_in} -> {x_expected, y_expected}
	std::vector<std::pair<std::pair<double,double>, std::pair<double,double>>> cases;
	
	auto project_to_circle = [&](double x_in, double y_in) {
	    double dx = x_in - x_c;
	    double dy = y_in - y_c;
	    double r2 = dx*dx + dy*dy;
	    if (r2 <= R_max*R_max) {
	        return std::make_pair(x_in, y_in); // inside, unchanged
	    } else {
	        double theta = std::atan2(dy, dx);
	        double x_out = x_c + R_max*std::cos(theta);
	        double y_out = y_c + R_max*std::sin(theta);
	        return std::make_pair(x_out, y_out);
	    }
	};
	
	// Define input points
	std::vector<std::pair<double,double>> inputs = {
	    {x_c + 0.5*R_max, y_c},           // inside
	    {x_c + R_max, y_c},               // on boundary
	    {x_c + 2*R_max, y_c},             // outside +x
	    {x_c - 2*R_max, y_c},             // outside -x
	    {x_c, y_c + 2*R_max},             // outside +y
	    {x_c, y_c - 2*R_max},             // outside -y
	    {x_c + 2*R_max, y_c + 2*R_max},   // outside 1st quadrant
	    {x_c - 2*R_max, y_c + 2*R_max},   // outside 2nd quadrant
	    {x_c - 2*R_max, y_c - 2*R_max},   // outside 3rd quadrant
	    {x_c + 2*R_max, y_c - 2*R_max}    // outside 4th quadrant
	};
	
	// Generate {input -> expected} pairs
	for (auto &pt : inputs) {
	    cases.push_back({pt, project_to_circle(pt.first, pt.second)});
	}

    for (size_t i = 0; i < cases.size(); ++i) 
	{
        double x = cases[i].first.first;
        double y = cases[i].first.second;
        wall.apply_boundary_condition(x, y);

        double x_expected = cases[i].second.first;
        double y_expected = cases[i].second.second;

        if (!(float_equality(x, x_expected, tol) && float_equality(y, y_expected, tol))) 
		{
            std::cerr << "Test " << i+1 << " failed: "
                      << "got (" << x << "," << y << "), "
                      << "expected (" << x_expected << "," << y_expected << ")\n";
            return false;
        }
    }

    return true;
}

