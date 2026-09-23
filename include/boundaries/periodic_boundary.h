#ifndef PERIODIC_BOUNDARY_H
#define PERIODIC_BOUNDARY_H

#include "boundary.h"
#include "../common.h"

/***************************************************** 
 *
 * class: Periodic_boundary
 *
 * Class for enforcing periodic boundary conditions
 *
 * These come in pairs
 * 
 ******************************************************/

class Periodic_boundary : public Boundary 
{
public:

	/// Set PBC in the direction wiith length _L and limits lim
	/// @param _L - distance between pairs of boundaries
	/// @param lim - real box limits (not necessarily 0 to _L)
	/// Sets everything to 0 to L for the ease of implementation
	Periodic_boundary(const double _L, std::vector<double> lim) : 
		L(_L), pos_0(lim.at(0)), pos_f(lim.at(1)), dL(pos_0) { }	

	/// Function for imposing BC on this boundary
	void apply_boundary_condition(double& x) override;

	/// Function for imposing BC on this boundary
	void apply_boundary_condition(double&, double&) override {;}

	/// Function for imposing particle-specific BC on this boundary
	void apply_boundary_condition_for_particle(double&, double) override {;}


protected:

	// Length in the direction periodic boundaries are imposed
	double L = 0.0;
	// Lower and upper limits - positions of the boundary pair
	double pos_0 = 0.0;
	double pos_f = 0.0;
	// Shift to have the domain go from 0 to L
	double dL = 0.0;
};

 #endif
