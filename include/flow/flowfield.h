#ifndef FLOWFIELD_H
#define FLOWFIELD_H

#include "flowimp.h"

/***************************************************** 
 * class: FlowField 
 *
 * Fetches and stores the fluid flow information as 
 *	a vector field constructed from discrete data  
 * 
 ******************************************************/

class FlowField : public FlowImp {
public:

	//
	// Constructors
	//

	// Default is disabled
	FlowField() = delete;	

	/// Creates a zero velocity field with Nx x Ny  nodes
	FlowField(const size_t Nx, const size_t Ny) 
		{ vel_x.resize(Nx*Ny, 0.0); vel_y.resize(Nx*Ny, 0.0); }

	//
	// Retrieving velocity values
	//

	/// x velocity component at position (x,y)
	double velocity_x(const double x, const double y)
		{ return 0.0; } // Placeholder
	/// y velocity component at position (x,y)
	double velocity_y(const double x, const double y)
		{ return 0.0; } // Placeholder

private:

	// Velocity field, x and y components
	std::vector<double> vel_x;
	std::vector<double> vel_y;

};

#endif
