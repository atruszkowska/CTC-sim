#ifndef CIRCLE_H
#define CIRCLE_H

#include "../../common.h"

/*************************************************************** 
 * class: Circle
 *
 * Defines a circular enclosure - well
 *
 ***************************************************************/

class Circle {
public:

	//
	// Constructors 	
	// 


 	/// Creates an object with default properties
	Circle() = default;

	//
	// Setup
	//

	/** 
	 * \brief Define dimensions 
	 * @param R - domain radius 
	 */
	inline void set_radius(const double R) 
	{ _R = R; }

	/** 
	 * \brief Define the center of the box 
	 * @param x_c - domain center, x (horizontal) coordinate
	 * @param y_c - domain center, y (vertical) coordinate
	 */
	inline void set_center(const double x_c, const double y_c) 
	{ _x_c = x_c; _y_c = y_c; }

	//	
	// Getters
	//
	
	/// Get the radius of the domain 
	inline double R() const { return _R; }
	/// Returns the coordinates of the center of the enclosure [x, y] 
    inline std::vector<double> get_center_coordinates() const { return {_x_c, _y_c}; }
	
private:

	// Radius of the domain (in units of length)
	double _R = 0.0;
	// Domain center
	double _x_c = 0.0, _y_c = 0.0;

};

#endif

