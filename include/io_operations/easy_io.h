#ifndef EASY_IO_H
#define EASY_IO_H

#include "file_handler.h"
#include "../common.h"

/*************************************************************** 
 * class: Easy_IO
 *
 * Interface for read and write operations 
 *
 * Manages IO operations for the EasyFlow simulations
 **************************************************************/

class Easy_IO
{
public:

	//
	// Constructors
	//
	
	/**
	 * \brief Default constructor - sets file names and writing mode to single file
	 * \details File will be located in the current working directory
	 */ 
	Easy_IO() = default;

	/** 
	 * \brief Creates an Easy_IO object
	 * @param name - name of the file or begining of the file names for file sequences
	 */
	Easy_IO(const std::string name) : fname(name) { } 

	//
	// Reading functionality 
	//
	
	/**
	 * \brief Read whitespace separated data from a file
	 * \details Data in a line needs to be separated by whitespace (tab, blank, ...)
	 * @returns std::vector of vectors of type T with the data in each vector being
	 * 		one line of the file with elements as extracted baed on whitespace
	 */
	template<typename T>
	std::vector<std::vector<T>> read_vector() const;

	//
	// Writing functionality 
	//

	/**
	 * \brief Write object information to a file
	 * \details Object information as defined by the overloaded << 
	 * operator. Truncates if the file exists.
	 * 
	 * @param object - instance of object of type T
	 */
	template<typename T>
	void write_object(const T& object) const;

	/**
	 * \brief Write a 1D vector to a file, each entry in a new line
	 * \details Truncates if the file exists. 
	 *
 	 * @param data - vector of type T
	 */
	template<typename T>
	void write_vector(const std::vector<T>& data) const;

	/**
	 * \brief Write a 2D vector to a file with delimiter delim
	 * \details Truncates if the file exists. 
     *
	 * @param data - nested vector of type T with dimensions dims[0] x dims[1] 
	 */
	template<typename T>
	void write_vector(const std::vector<std::vector<T>>& data) const;
	
	//
	// Destructor
	//
	
	~Easy_IO() = default;

private:
	std::string fname = {"dflt_easy_io_file.txt"};
};

//
// Reading functionality 
//

template<typename T>
std::vector<std::vector<T>> Easy_IO::read_vector() const
{
	File_handler file(fname);
	std::fstream &in = file.get_stream();
	std::string line;
   	T elem;
	std::vector<std::vector<T>> output;

	while (std::getline(in, line)){
		std::istringstream data_row(line);
		std::vector<T> temp;
		while (data_row >> elem) {
			temp.push_back(elem);
		}
		output.push_back(temp);
	}
	return output;
}

//
// Writing functionality 
//

template<typename T>
void Easy_IO::write_object(const T& object) const
{
	File_handler file(fname, std::ios_base::out | std::ios_base::trunc);
	std::fstream &fout = file.get_stream();
	fout << object << '\n';
}

template<typename T>
void Easy_IO::write_vector(const std::vector<T>& data) const
{
	File_handler file(fname, std::ios_base::out | std::ios_base::trunc);
	std::fstream &fout = file.get_stream();
	for (const auto& entry : data) {
		fout << entry << '\n';		
	}
}

template<typename T>
void Easy_IO::write_vector(const std::vector<std::vector<T>>& data) const
{
	File_handler file(fname, std::ios_base::out | std::ios_base::trunc);
	std::fstream &out = file.get_stream();
	size_t nrows = data.size();
	for (size_t i = 0; i<nrows; i++){
		std::copy(data.at(i).begin(), data.at(i).end(),	std::ostream_iterator<T>(out, " "));
		out << '\n';
	}	
}		

#endif
