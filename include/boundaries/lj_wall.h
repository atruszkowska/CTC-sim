#ifndef LJ_WALL_H
#define LJ_WALL_H

#include "boundary.h"
#include "../common.h"

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

class Lj_wall : public Boundary 
{
public:

	/// Constructor 
	// @lim  – real box limits in the wall direction (pos0, posF)
	// @_sigma   – LJ length scale (particle radius)     (m)
	// @_epsilon – well depth (energy units)
	// @_alpha - factor that multiplies epsilon (unitless)
	// @_m - first power
	// @_n - second power
	// @_dt - time step
	// @_m_cell - mass of the cell
	// @_Fc - the force that acts on the second coordinate - nominal magnitude and direction
	Lj_wall(const double _L, const std::vector<double> lim, 
		const double _sigma, const double _epsilon, 
		const double _alpha, const double _Fc, const int _m, const int _n, 
		const double _dt, const double _m_cell) : 
		L(_L), pos_0(lim.at(0)), pos_f(lim.at(1)), 
		dL(pos_0), sigma(_sigma), epsilon(_epsilon), alpha(_alpha), Fc(_Fc),
		m(_m), n(_n), dt(_dt), m_cell(_m_cell){}
	
	// Function for imposing BC on this boundary
	void apply_boundary_condition(double& x) override;	

	// Function for imposing BC on this boundary - both coordinates
	void apply_boundary_condition(double& x, double& y) override;

	void apply_boundary_condition_for_particle(double&, double) override {;}

	// Generate the entire LJ force over a predefined range
	// @x0 - initial position
	// @xf - final position
	// n_samples - number of values to generate (evenly spaced)
	const std::vector<double> compute_lj_force(const double x0, const double xf, const int n_samples);	

protected:

	// Length in the direction LJ wall boundaries are imposed
	double L = 0.0;
	// Lower and upper limits - positions of the boundary pair
	double pos_0 = 0.0;
	double pos_f = 0.0;
	// Shift to have the domain go from 0 to L
	double dL = 0.0;
	
	// LJ parameters
    double sigma =0.0;
    double epsilon = 0.0;
	double alpha = 0.0;
	double Fc = 0.0;
    int m = 0.0;
	int n = 0.0;

	// Other parameters
	double dt = 0.0;
	double m_cell = 0.0;

};

#endif
