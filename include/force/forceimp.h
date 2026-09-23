#ifndef FORCEIMP_H
#define FORCEIMP_H

#include "../common.h"

/***************************************************** 
 * class: ForceImp 
 *
 * Abstract base class for external force field 
 *	implementation  
 * 
 ******************************************************/

class ForceImp {
public:
	// Functions that return forces that act on 
	// x and y position components
	virtual double force_on_x(double) = 0; 	
	virtual double force_on_y(double) = 0;

	// Destructors
	virtual ~ForceImp() = default;
}; 

#endif
