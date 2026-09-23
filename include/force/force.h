#ifndef FORCE_H
#define FORCE_H

#include "../common.h"
#include <memory>

#include "forcefunc.h"

class ForceFunction;

/***************************************************** 
 * class: Force
 *
 * Stores, manages, and provides the information about
 * external force field.
 * 
 ******************************************************/

class Force {
public:

 	// Function that returns an external force value 
	// at given position x or y
	typedef double(*eforce)(double);

	//
	// Constructors
	//

	/// No velocity field at all
	Force() = default;

	/// Spatial relation for force components
	/// @param F_x - function that takes x coordinate
	///					and returns the force that acts on x
	/// @param F_y - function that takes y coordinate
	///					and returns the force that acts on y
	Force(eforce F_x, eforce F_y) : 
		imp_ptr(new ForceFunc(F_x, F_y)) { }

	/// Copy constructor - disabled 
	Force(const Force&) = delete;

	/// Assignment operator - disabled 
	const Force& operator= (const Force&) = delete;

	//
	// Retrieving velocity values
	//

	/// Force acting on the x component of position (x,y)
	double force_on_x(const double x)
		{ return imp_ptr->force_on_x(x); } 
	/// Force acting on the y component of position (x,y)
	double force_on_y(const double y)
		{ return imp_ptr->force_on_y(y); }

private:

	// Pointer to class with velocity field implementation
	std::unique_ptr<ForceImp> imp_ptr{nullptr};

};

#endif
