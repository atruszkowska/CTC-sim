#include "../../include/distances/periodic_distance.h"
#include "../../include/distances/cartesian_distance.h"
#include "../../include/test_utils.h"

/***************************************************** 
 *
 * Test suite for the Distance class hierarchy 
 *
 *****************************************************/

// Tests
bool periodic_distance_test();
bool cartesian_distance_test();

int main()
{
	test_pass(periodic_distance_test(), "Periodic distance");
	test_pass(cartesian_distance_test(), "Cartesian distance");
}

bool periodic_distance_test()
{
	double tol = 1e-5;
	std::vector<double> x1{2.5, 0.0, -1.3, 100.9};
	std::vector<double> x2{1.5, 1.0, -10.8, 100.9};
	std::vector<std::vector<double>> exp_dist = {{1.0, 1.0, -7.5, 0.0}, {1.0, 1.0, 0.5, 0.0}};
	std::vector<double> L{2.0, 10.0};
	double dist = 0.0;

	for (size_t ib = 0; ib < L.size(); ++ib)
	{
		Periodic_distance periodic(L.at(ib));
	
		for (size_t ij = 0; ij < x1.size(); ++ij)
		{
			dist = periodic.compute_distance(x1.at(ij), x2.at(ij));
			if (!(float_equality(dist, exp_dist.at(ib).at(ij), tol))) 
			{
				std::cerr << "Computed distance " << dist  
						  << " does not match expected " 
						  << exp_dist.at(ib).at(ij) << std::endl;
				return false;
			} 	
		}

	}
	return true;
}

bool cartesian_distance_test()
{
	double tol = 1e-5;
	std::vector<double> x1{2.5, 0.0, -1.3, 100.9};
	std::vector<double> x2{1.5, 1.0, -10.8, 100.9};
	std::vector<double> exp_dist{1.0, -1.0, 9.5, 0.0};
	double dist = 0.0;	

	Cartesian_distance cartesian;

	for (size_t ij = 0; ij < x1.size(); ++ij)
	{
		dist = cartesian.compute_distance(x1.at(ij), x2.at(ij));
		if (!(float_equality(dist, exp_dist.at(ij), tol))) 
		{
			std::cerr << "Computed distance " << dist  
					  << " does not match expected " 
					  << exp_dist.at(ij) << std::endl;
			return false;
		} 	
	}

	return true;
}
