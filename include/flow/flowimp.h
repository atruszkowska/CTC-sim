#ifndef FLOWIMP_H
#define FLOWIMP_H

#include "../common.h"

/***************************************************** 
 * class: FlowImp 
 *
 * Abstract base class for external flow field 
 *	implementation  
 * 
 ******************************************************/

class FlowImp {
public:
	// Functions that return x and y velocity components
	virtual double velocity_x(double, double) = 0; 	
	virtual double velocity_y(double, double) = 0;

	// Destructors
	virtual ~FlowImp() = default;
}; 

#endif
