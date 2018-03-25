#include "tools.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

void fprint_vector(std::string filename, const std::vector<double> &field)
{
	std::ofstream output_file;
	output_file.open(filename, std::ios::out);
	
	for (auto data : field)
	{
		output_file << data << std::endl;
	}
	
	output_file.close();
	
	return;
}

void fprint_vector(const std::string prefix, int iter, const std::vector<double> &field)
{
	std::ostringstream stringStream;
	stringStream << prefix << '_' << std::setfill('0') << std::setw(4) << iter << ".txt";
	
	std::ofstream output_file;
	output_file.open(stringStream.str(), std::ios::out);
	
	for (auto data : field)
	{
		output_file << data << std::endl;
	}
	
	output_file.close();
	
	return;
}