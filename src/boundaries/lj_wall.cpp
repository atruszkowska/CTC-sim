#include "../../include/boundaries/lj_wall.h"

/***************************************************** 
 *
 * class: LJ_boundary
 *
 * Class for enforcing Lennard-Jones potential as a 
 * wall boundary condition
 *
 * Implements a pair of confining walls each exerting
 * a shifted M-N LJ potential on particle centres:
 * 
 * The potential:
 *
 * U(r) = alpha*epsilon*(sigma^M/r^M - sigma^N/r^N)
 *
 * The resulting force:
 *
 * F(r) = alpha*epsilon(M*sigma^M/r^(M+1) - N*sigma^N/r^(N+1))
 *
 * Where r is the distance between particle center and 
 * the wall. 
 *
 ******************************************************/

// Function for imposing BC on this boundary
void Lj_wall::apply_boundary_condition(double& x)
{
	// Compute distance with each wall
	const double dist_0 = std::fabs(x - pos_0);
	const double dist_f = std::fabs(x - pos_f);

	// Compute forces with each wall and net force
	// F0 sign is flipped to achieve the desired displacement 
	const double F0 = alpha*epsilon*(m*std::pow(sigma, m)/std::pow(dist_0, (m+1)) 
						- n*std::pow(sigma, n)/std::pow(dist_0, (n+1)));
	const double Ff = -alpha*epsilon*(m*std::pow(sigma, m)/std::pow(dist_f, (m+1)) 
						- n*std::pow(sigma, n)/std::pow(dist_f, (n+1)));
	const double F_net = F0 + Ff;

	// Compute the displacement from the net force, correct the position  
	x += (F_net*dt*dt/m_cell);
}

// Function for imposing BC on this boundary
void Lj_wall::apply_boundary_condition(double& x, double& y)
{
	// Compute distance with each wall
	const double dist_0 = std::fabs(x - pos_0);
	const double dist_f = std::fabs(x - pos_f);

	// Compute forces with each wall and net force
	// F0 sign is flipped to achieve the desired displacement 
	const double F0 = alpha*epsilon*(m*std::pow(sigma, m)/std::pow(dist_0, (m+1)) 
						- n*std::pow(sigma, n)/std::pow(dist_0, (n+1)));
	const double Ff = -alpha*epsilon*(m*std::pow(sigma, m)/std::pow(dist_f, (m+1)) 
						- n*std::pow(sigma, n)/std::pow(dist_f, (n+1)));
	const double F_net = F0 + Ff;

	// Compute the displacement from the net force, correct the position  
	// This is the coordinate aligned with the of the LJ wall
	x += (F_net*dt*dt/m_cell);

	// This is the coordinate parallel to the LJ wall (and its normal)
	y += (F_net*dt*dt/m_cell)*Fc;

}

// Generate and save the entire LJ force over a predefined range
const std::vector<double> Lj_wall::compute_lj_force(const double x0, const double xf, const int n_samples)
{
	const double dx = (xf - x0)/n_samples;
	std::vector<double> lj_force = {};
	double xi = x0; 

	while (xi <= xf)
	{
		lj_force.push_back(alpha*epsilon*(m*std::pow(sigma, m)/std::pow(xi, (m+1)) 
						- n*std::pow(sigma, n)/std::pow(xi, (n+1))));
		xi += dx; 
	}

	return lj_force;
}
