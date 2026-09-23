#include "../../../include/flow/flow.h"
#include "../../../include/utils.h"
#include "../../../include/test_utils.h"

// Tests
bool simple_profile();

int main()
{
	test_pass(simple_profile(), "Basic xy function");
}

bool simple_profile()
{
	// For comparison
	double tol = 1e-5;
	// x component
	auto vel_x = [](double x, double y) { return 2*x + y; }; 
	// y component
	auto vel_y = [](double x, double y) { return 2*x + 3*y; };
	// Values to test with
	std::vector<std::vector<double>> coords{{0,0}, {-10, 1}, {10, -1}, 
										{1.7, 2.8}, {-3, -5.1}};

	// Implementation class 
	Flow flow(vel_x, vel_y);

	// Test
	for (size_t iv = 0; iv < coords.size(); ++iv)
	{
		if (!(float_equality(
				vel_x(coords.at(iv).at(0), coords.at(iv).at(1)),
				flow.velocity_x(coords.at(iv).at(0), coords.at(iv).at(1)), tol)))
		{
			std::cerr << "x velocity " 
				<< flow.velocity_x(coords.at(iv).at(0), coords.at(iv).at(1))  
				<< " does not match expected " 
				<< vel_x(coords.at(iv).at(0), coords.at(iv).at(1)) << std::endl;
			return false;
		}
		if (!(float_equality(
				vel_y(coords.at(iv).at(0), coords.at(iv).at(1)),
				flow.velocity_y(coords.at(iv).at(0), coords.at(iv).at(1)), tol)))
		{
			std::cerr << "y velocity " 
				<< flow.velocity_y(coords.at(iv).at(0), coords.at(iv).at(1))  
				<< " does not match expected " 
				<< vel_y(coords.at(iv).at(0), coords.at(iv).at(1)) << std::endl;
			return false;
		}
	}

	return true;
}
