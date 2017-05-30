#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "TVD_scheme.h"

void print_field(std::string filename, std::vector<double> &field);

void get_field_time(std::vector<std::vector<double>> &timelayer_fields, std::vector<double> &field, int time_layers_cnt, int time_layer);
void add_field_time(std::vector<std::vector<double>> &timelayer_fields, std::vector<double> &field, int time_layers_cnt, int time_layer);

void make_GL_coeff_vec(std::vector<double> &GL_vec, double alpha, int n);

int main()
{
	const int Nx = 100;
	const double dx = 1.;
	
	double dt = 0.01;
	double finish_time = 10.;
	
	// frac
	double t_len = 5.;
	int    time_layers_cnt = int(t_len / dt);
	std::vector<std::vector<double>> timelayer_fields(time_layers_cnt, std::vector<double>(Nx, 0.));
	
	// fields
	std::vector<double> vel(Nx - 1, 2.);
	std::vector<double> field(Nx, 1.);
	std::vector<double> field_new(Nx, 0.);
	
	// initial field values
	const int ix_start = 10, ix_end = 20;
	for(int ix = ix_start; ix < ix_end; ix++)
	{
		field[ix] = 5.;
	}
	timelayer_fields[0] = field;
	
	// solver
	TVD_scheme solve_scheme(Nx, dx);
	
	// time counters
	double time = 0.;
	int it = 0;
	
	const double alpha = 0.7;
	std::vector<double> GL_coeffs;
	make_GL_coeff_vec(GL_coeffs, alpha, time_layers_cnt);
	print_field(std::string("GL_coeffs.txt"), GL_coeffs);
	
	// time cycle
	while(time < finish_time)
	{
		// solve
		solve_scheme.solve_transfer_explicitly(vel, field, field_new, dt);
		
		// time counters ++
		field = field_new;
		time += dt;
		it++;
		std::cout << "Calculated " << time << " secs" << std::endl;
		
		// save last field
		add_field_time(timelayer_fields, field_new, time_layers_cnt, it);
		
		// file print
		std::ostringstream stringStream;
		stringStream << "field_" << std::setfill('0') << std::setw(int(log10(finish_time / dt)) + 1) << it << ".txt";
		print_field(stringStream.str(), field_new);
	}
	
	// print final result
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

void get_field_time(std::vector<std::vector<double>> &timelayer_fields, std::vector<double> &field, int time_layers_cnt, int time_layer)
{
	field = timelayer_fields[time_layer % time_layers_cnt];
	
	return;
}

void add_field_time(std::vector<std::vector<double>> &timelayer_fields, std::vector<double> &field, int time_layers_cnt, int time_layer)
{
	timelayer_fields[time_layer % time_layers_cnt] = field;
	
	return;
}

void make_GL_coeff_vec(std::vector<double> &GL_vec, double alpha, int n)
{
	GL_vec.resize(n + 1, 0.);
	
	GL_vec[0] = 1.;
	for(int i = 1; i < n + 1; i++)
	{
		double mult = (alpha - i + 1) / i;
		GL_vec[i] = GL_vec[i - 1] * mult;
	}
	
	for(int i = 1; i < n + 1; i+=2)
		GL_vec[i] *= -1.;
	
	return;
}
