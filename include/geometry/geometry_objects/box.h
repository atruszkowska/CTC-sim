#ifndef BOX_H
#define BOX_H

#include "../../common.h"

/*************************************************************** 
 * class: Box 
 *
 * Defines a rectangular enclosure - box
 *
 ***************************************************************/

class Box {
public:

	//
	// Constructors 	
	// 


 	/// Creates a Box object with default properties
	Box() = default;

	//
	// Setup
	//

	/** 
	 * \brief Define box dimensions 
	 * @param Lx - domain dimension in x (horizontal coordinate) 
	 * @param Ly - domain dimension in y (vertical coordinate)
	 */
	inline void set_dimensions(const double Lx, const double Ly) 
	{ _Lx = Lx; _Ly = Ly; }

	/** 
	 * \brief Define the center of the box 
	 * @param x_c - domain center, x (horizontal) coordinate
	 * @param y_c - domain center, y (vertical) coordinate
	 */
	inline void set_center(const double x_c, const double y_c) 
	{ _x_c = x_c; _y_c = y_c; }

	/// Find and store object bounds
	inline void compute_bounds()
		{ x_0 = _x_c - _Lx/2; y_0 = _y_c - _Ly/2; 
	  		x_f = _x_c + _Lx/2; y_f = _y_c + _Ly/2; }

	//	
	// Getters
	//
	
	/// Get the length of the domain in x (horizontal) direction
	inline double Lx() const { return _Lx; }
	/// Get the length of the domain in y (vertical) direction
    inline double Ly() const { return _Ly; }
	/// Returns the lower and upper bounds of a box enclosure [x0, y0, xf, yf] 
    inline std::vector<double> get_limits() const { return {x_0, y_0, x_f, y_f}; }
	
private:

	// Dimensions of the domain (in units of length)
	double _Lx = 0.0, _Ly = 0.0;
	// Domain center
	double _x_c = 0.0, _y_c = 0.0;
	// Domain bounds
	double x_0 = 0.0, y_0 = 0.0, x_f = 0.0, y_f = 0.0;

};

#endif

