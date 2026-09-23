#ifndef WALL_BOUNDARY_H
#define WALL_BOUNDARY_H

#include "boundary.h"
#include "../common.h"

/***************************************************** 
 *
 * class: Wall_boundary
 *
 * Class for enforcing wall boundary conditions
 *
 * This also acts in pairs - it assumes each wall has
 * another wall opposite to it (to be generalized)
 *
 ******************************************************/

class Wall_boundary : public Boundary 
{
public:

	/// A pair of walls in a direction
	/// Sets everything to 0 to L for the ease of implementation
	/// @param _L - distance between pairs of boundaries
	/// @param lim - real box limits (not necessarily 0 to _L)
	/// @param _R_min - minimum allowable distance with the wall
	Wall_boundary(const double _L, std::vector<double> lim, 
		const double _R_min) : 
		L(_L), pos_0(lim.at(0)), pos_f(lim.at(1)), 
		dL(pos_0), R_min(_R_min) { }
	
	Wall_boundary(const double _L, std::vector<double> lim) : 
		L(_L), pos_0(lim.at(0)), pos_f(lim.at(1)), 
		dL(pos_0) { }

	/// Function for imposing BC on this boundary
	void apply_boundary_condition(double& x) override;

	/// Function for imposing BC on this boundary
	void apply_boundary_condition(double&, double&) override {;}

	/// Function for imposing particle-specific BC on this boundary
	void apply_boundary_condition_for_particle(double& x, double rmin) override;

protected:

	// Length in the direction wall boundaries are imposed
	double L = 0.0;
	// Lower and upper limits - positions of the boundary pair
	double pos_0 = 0.0;
	double pos_f = 0.0;
	// Shift to have the domain go from 0 to L
	double dL = 0.0;
	// Minimum allowable distance with the wall
	double R_min = 0.0;

};

#endif
