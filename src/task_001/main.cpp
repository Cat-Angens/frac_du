#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "TVD_scheme.h"

void print_field(std::string filename, std::vector<double> &field);


int main()
{
	const int Nx = 100;
	const double dx = 1.;
	
	double dt = 0.01;
	double finish_time = 10.;
	
	std::vector<double> vel(Nx - 1, 2.);
	std::vector<double> field(Nx, 1.);
	std::vector<double> field_new(Nx, 0.);
	
	const int ix_start = 10, ix_end = 20;
	for(int ix = ix_start; ix < ix_end; ix++)
	{
		field[ix] = 5.;
	}
	
	TVD_scheme solve_scheme(Nx, dx);
	
	double time = 0.;
	int it = 0;
	while(time < finish_time)
	{
		solve_scheme.solve_transfer_explicitly(vel, field, field_new, dt);
		field = field_new;
		time += dt;
		it++;
		std::cout << "Calculated " << time << " secs" << std::endl;
		
		std::ostringstream stringStream;
		stringStream << "field_" << std::setfill('0') << std::setw(int(log10(finish_time / dt)) + 1) << it << ".txt";
		print_field(stringStream.str(), field_new);
	}
	
	std::string filename("out.txt");
	print_field(filename, field_new);
	
	return 0;
}

void print_field(std::string filename, std::vector<double> &field)
{
	std::ofstream output_file;
	output_file.open(filename, std::ios::out);
	
	for(auto data : field)
	{
		output_file << data << std::endl;
	}
	
	return;
}

