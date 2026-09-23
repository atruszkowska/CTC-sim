#ifndef FLOWFUNC_H
#define FLOWFUNC_H

#include "flowimp.h"

/***************************************************** 
 * class: FlowFunc 
 *
 * Fetches the fluid flow information from 
 *	a user-defined function  
 * 
 ******************************************************/

class FlowFunc : public FlowImp {
public:

	// Function that returns a velocity value 
	// at a given position x and y
	typedef double(*vel)(double, double);

	//
	// Constructors
	//

	// Default is disabled
	FlowFunc() = delete;	

	/// Stores user-defined velocity relations
	FlowFunc(vel v_x, vel v_y) : vel_x(v_x), vel_y(v_y) { }

	//
	// Retrieving velocity values
	//

	/// x velocity component at position (x,y)
	double velocity_x(const double x, const double y)
		{ return vel_x(x, y); } 
	/// y velocity component at position (x,y)
	double velocity_y(const double x, const double y)
		{ return vel_y(x, y); } 

private:

	// Velocity field, x and y components
	vel vel_x;
	vel vel_y;

};

#endif
