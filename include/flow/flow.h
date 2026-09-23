#ifndef FLOW_H
#define FLOW_H

#include "../common.h"
#include <memory>

#include "flowfield.h"
#include "flowfunc.h"

class FlowField;
class FlowFunction;

/***************************************************** 
 * class: Flow 
 *
 * Stores, manages, and provides the information about
 * external fluid flow.
 * 
 ******************************************************/

class Flow {
public:

	// Function that returns a velocity value 
	// at a given position x and y
	typedef double(*vel)(double, double);

	//
	// Constructors
	//

	/// No velocity field at all
	Flow() = default;

	/// No flow, zero velocity
	Flow(const size_t Nx, const size_t Ny) : 
		imp_ptr(new FlowField(Nx, Ny)) {};

	/// Spatial relation for velocity components
	/// @param v_x - function that takes x and y spatial coordinates
	///					and returns x velocity component
	/// @param v_y - function that takes x and y spatial coordinates
	///					and returns y velocity component
	Flow(vel v_x, vel v_y) : 
		imp_ptr(new FlowFunc(v_x, v_y)) { }

	/// Copy constructor - disabled 
	Flow(const Flow&) = delete;

	/// Assignment operator - disabled 
	const Flow& operator= (const Flow&) = delete;

	//
	// Retrieving velocity values
	//

	/// x velocity component at position (x,y)
	double velocity_x(const double x, const double y)
		{ return imp_ptr->velocity_x(x, y); }
	/// y velocity component at position (x,y)
	double velocity_y(const double x, const double y)
		{ return imp_ptr->velocity_y(x, y); }

private:

	// Pointer to class with velocity field implementation
	std::unique_ptr<FlowImp> imp_ptr{nullptr};

};

#endif
