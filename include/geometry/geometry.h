#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "../common.h"
#include "geometry_include.h"

/*************************************************************** 
 * class: Geometry
 *
 * Creation, storage, and access of model geometry
 *
 ***************************************************************/

class Geometry {
public:

	//
	// Constructors 	
	// 


 	/// Creates a Geometry object with default properties
	Geometry() = default;

	//
	// Domains/enclosures 
	//

	/** 
	 * \brief Create a rectangular domain
	 * @param Lx - domain dimension in x (horizontal coordinate) 
	 * @param Ly - domain dimension in y (vertical coordinate)
	 * @param x_c - domain center, x (horizontal) coordinate
	 * @param y_c - domain center, y (vertical) coordinate
	 */
	inline void make_box(const double Lx, const double Ly, 
			const double x_c, const double y_c) 
	{ enclosure.set_dimensions(Lx, Ly); 
	  	  enclosure.set_center(x_c, y_c); 
		  enclosure.compute_bounds(); }

	/** 
	 * \brief Create a circular domain
	 * @param R - domain radius 
	 * @param x_c - domain center, x (horizontal) coordinate
	 * @param y_c - domain center, y (vertical) coordinate
	 */
	inline void make_circle(const double R, const double x_c, const double y_c) 
		{ well.set_radius(R); well.set_center(x_c, y_c); }

	//	
	// Getters
	//
	
	/// Returns the length of the domain in x (horizontal) direction
	inline double get_x_length() const { return enclosure.Lx(); }
	/// Returns the length of the domain in y (vertical) direction
    inline double get_y_length() const { return enclosure.Ly(); }
	/// Returns the lower and upper limits of a box enclosure [x0, y0, xf, yf]
    inline std::vector<double> get_box_limits() const { return enclosure.get_limits(); }
	/// Returns the radius of a circular enclosure
	inline double get_radius() const { return well.R(); }
	/// Returns the center coordinates of a circular enclosure
	inline std::vector<double> get_center() const { return well.get_center_coordinates(); }

private:

	// Domain - rectangular enclosure 
	Box enclosure;
	// Domain - circular enclosure
	Circle well;

};

#endif

