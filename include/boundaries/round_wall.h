#ifndef ROUND_WALL_H
#define ROUND_WALL_H

#include "boundary_xy.h"
#include "../common.h"

/***************************************************** 
 *
 * class: Round_wall
 *
 * Class for enforcing wall boundary conditions inside
 * of a circular encolosure
 *
 * 
 ******************************************************/

class Round_wall : public Boundary_xy 
{

public:

	/// Defines the enclosure 	
	/// @param _R - enclosure radius
	/// @param coords - real coordinates of the center of the enclosure, (x, y)
	/// @param _R_min - minimum allowable distance with the wall
	Round_wall(const double _R, std::vector<double> coords, 
		const double _R_min) : 
		R(_R), x_c(coords.at(0)), y_c(coords.at(1)), 
		R_max(R - _R_min) { }
	
	Round_wall(const double _R, std::vector<double> coords) : 
		R(_R), x_c(coords.at(0)), y_c(coords.at(1)), 
		R_max(R) { }

	/// Function for imposing BC on this boundary
	/// The arguments are x and y coordinates of a particle
	void apply_boundary_condition(double&, double&) override;
	/// Function for imposing particle-specific BC on this boundary
	void apply_boundary_condition_for_particle(double&, double&, double) override;

	
protected:

	// Enclosure radius 
	double R = 0.0;
	// Real coordinates of the center of the enclosure
	double x_c = 0.0;
	double y_c = 0.0;
	// Maximum allowable radial position 
	double R_max = 0.0;

};

#endif
