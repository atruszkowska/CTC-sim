#ifndef LOAD_PARAMETERS_H
#define LOAD_PARAMETERS_H

#include "file_handler.h"
#include "../common.h"


/***********************************************************************************************
 *
 * class: LoadParameters
 *
 * Class for reading and storing parameters from file
 *
 *
 * ****************************************************************************************/

class LoadParameters
{
public:

	LoadParameters() = default;


	/**
	 *  \brief Read parameters from file, store them as a map
	 *  \details In the file parameter name can be multiple words
	 *  		but it needs to be in a seperate line with a // comment
	 *  		type as the first element of the line; it is followed
	 *  		by a numeric value in the next line.
	 *
	 *
	 *  @param infile - input file with parameters
	 */
	std::map<std::string, double> load_parameter_map(const std:: string infile);

};

#endif
