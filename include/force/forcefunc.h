#ifndef FORCEFUNC_H
#define FORCEFUNC_H

#include "forceimp.h"

/***************************************************** 
 * class: ForceFunc 
 *
 * Fetches the force field information from 
 *	a user-defined function  
 * 
 ******************************************************/

class ForceFunc : public ForceImp {
public:

 	// Function that returns an external force value 
	// at given position x or y
	typedef double(*eforce)(double);

	//
	// Constructors
	//

	// Default is disabled
	ForceFunc() = delete;	

	/// Stores user-defined force relations
	ForceFunc(eforce _F_x, eforce _F_y) : F_x(_F_x), F_y(_F_y) { }

	//
	// Retrieving force values
	//

	/// Force acting on the x component of position (x,y)
	double force_on_x(const double x)
		{ return F_x(x); } 
	/// Force acting on the y component of position (x,y)
	double force_on_y(const double y)
		{ return F_y(y); } 

private:

	// Force field, x and y components
	eforce F_x;
	eforce F_y;

};

#endif
