#include "../../../include/force/forcefunc.h"
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
	auto force_x = [](double x) { return 2*x; }; 
	// y component
	auto force_y = [](double y) { return 3*y + 1; };
	// Values to test with
	std::vector<std::vector<double>> coords{{0,0}, {-10, 1}, {10, -1}, 
										{1.7, 2.8}, {-3, -5.1}};

	// Implementation class 
	ForceFunc force(force_x, force_y);

	// Test
	for (size_t iv = 0; iv < coords.size(); ++iv)
	{
		if (!(float_equality(
				force_x(coords.at(iv).at(0)),
				force.force_on_x(coords.at(iv).at(0)), tol)))
		{
			std::cerr << "x force " 
				<< force.force_on_x(coords.at(iv).at(0))  
				<< " does not match expected " 
				<< force_x(coords.at(iv).at(0)) << std::endl;
			return false;
		}
		if (!(float_equality(
				force_y(coords.at(iv).at(1)),
				force.force_on_y(coords.at(iv).at(1)), tol)))
		{
			std::cerr << "y force " 
				<< force.force_on_y(coords.at(iv).at(1))  
				<< " does not match expected " 
				<< force_y(coords.at(iv).at(1)) << std::endl;
			return false;
		}
	}

	return true;
}
