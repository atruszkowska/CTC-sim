#include "../../include/io_operations/load_parameters.h"


/***********************************************************************************************
 * class: LoadParameters
 *
 * Class for reading and storing parameters from file
 *
 * ********************************************************************************************/

// Read parameters from file, store them as a map
std::map<std::string, double> LoadParameters::load_parameter_map(const std::string infile)
{
	File_handler file(infile);
	std::fstream &in = file.get_stream();
	std::string line;

	std::map<std::string, double> parameters;

	std::string tag = {}, word = {};

	while (std::getline(in, line)){

		std::istringstream data_row(line);

		// If a comment with a tag, save all words
		if (line.find("//") != std::string::npos){
			// Skip the // entry and blank line,
			// then collect the rest
			data_row >> word;
			while (data_row >> word)
				tag += (word + " ");
			// Remove the extra blank character
			if (!tag.empty())
				tag.pop_back();
			continue;
		}

		// If not - collect the number, make a map entry,
		// and reset the tag
		data_row >> word;
		parameters[tag] = std::stod(word);
		tag.clear();

	}
	return parameters;
}
