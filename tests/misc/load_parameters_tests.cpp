#include "../../include/test_utils.h"
#include <string>
#include "../../include/io_operations/load_parameters.h"


/*****************************************************************************************
 * Suite for testing LoadParameters class
 
 ****************************************************************************************/

// Tests
bool read_parameters_test();

// Supporting functions
bool equal_maps(std::map<std::string, double>, std::map<std::string, double>);

int main() 
{
	test_pass(read_parameters_test(), "Load particle parameters");
}

// Test for loading particle parameters
bool read_parameters_test()
{
	std::map<std::string, double> expected = 
			{{"v", 0.03}, {"radius", 0.5}, {"eta", 0.3}, {"dmin", 0.15}};
	std::map<std::string, double> loaded = {};

	LoadParameters ldparam;
	loaded = ldparam.load_parameter_map("test_data/ldp_input.txt");

	return equal_maps(expected, loaded);
}

bool equal_maps(std::map<std::string, double> expected, std::map<std::string, double> loaded)
{
	if (expected.size() != loaded.size()) 
	{
		return false;
	}
	for (const auto& entry : loaded) {
		const auto iter = expected.find(entry.first);
		if (iter == expected.end()){
			return false;
		} else {
			if (!float_equality<double>(iter->second, entry.second, 1e-5)) 
				return false;
		}
	}

	return true;
}
