#include "../../include/io_operations/file_handler.h"

/*************************************************************** 
 * class: File_handler
 *
 * Creates and manages input and output file streams
 *
 * Opens and maintains a file stream, then closes it when
 * calling the destructor. 
 **************************************************************/

void File_handler::open_file()
{
	file.open(fname, mode);
	if (!file.is_open()){
		std::cerr << "Error opening file " << fname << std::endl;
		throw std::ios_base::failure(std::strerror(errno));
	}
}

File_handler::~File_handler()
{
	// This is all that is being checked since std::getline sets failbit on EOF
	// https://stackoverflow.com/a/7855777/2763915
	// Leaving it like this for now - maybe extend in the future
	if (file.bad())
		std::cerr << "Error reading or writing to file " << fname 
				  << "\nDetected badbit set i.e. broken stream." << std::endl;
	file.close();
}		

