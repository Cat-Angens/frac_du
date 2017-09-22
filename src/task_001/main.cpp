#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "TVD_scheme.h"
#include "mat_solve_Thomas.h"

double get_Gamma(const double x);

void print_field(std::string filename, std::vector<double> &field);

void get_field_time(std::vector<std::vector<double>> &timelayer_fields, std::vector<double> &field, int time_layers_cnt, int time_layer);
void add_field_time(std::vector<std::vector<double>> &timelayer_fields, std::vector<double> &field, int time_layers_cnt, int time_layer);

void make_GL_coeff_vec(std::vector<double> &GL_vec, double alpha, int n);

// f(x, t)
double get_source(const double x, const double time)
{
	return 0.;
}

// phi(x)
double get_initial_field(const double x)
{
	if (x > 10. && x < 20.)
		return 1. - (x - 15.) * (x - 15.) / 25.;
	return 0.;
}

// phi'(x)
double get_initial_field_derivative(const double x)
{
	if (x > 10. && x < 20.)
		return -2. * (x - 15.) / 25.;
	return 0.;
}

// psi(t)
double get_border_value(const double time)
{
	return 0.;
}

void get_initial_field_vector(const double dx, const int Nx, std::vector<double> &field);
void get_source_vector(const double time, const double dx, const int Nx, std::vector<double> &source);


int main()
{
	const int Nx = 250;
	const double dx = 0.4;
	
	double dt = 0.1;
	double finish_time = 10.;
	
	const double alpha = 0.99;
	const double gamma_alpha = get_Gamma(alpha);
	
	// memory effects variables
	double t_len = 5.;
	int    time_layers_cnt = std::max(2, int(t_len / dt)); // minimum 2 layers - current and previous
	std::vector<std::vector<double>> timelayer_fields(time_layers_cnt, std::vector<double>(Nx, 0.));
	std::cout << "Time memory " << t_len << " secs (" << time_layers_cnt << " time steps)" << std::endl;
	
	////test Thomas method
	//// matrix 5x5
	//std::vector<std::vector<double>> matrix(3, std::vector<double>());
	//matrix[0] = {0., -0.1, -0.1, -0.1, -0.1};
	//matrix[1] = {1., 1., 1., 1., 1.};
	//matrix[2] = {0.1, 0.1, 0.1, 0.1, 0.};
	//std::vector<double> right_part = {1.2, 2.2, 3.2, 4.2, 4.6};
	//std::vector<double> solution(5, 0.);
	//mat_solve_Thomas solver(5);
	//solver.solve(matrix, right_part, solution);
	//
	//std::cout << "Solution:\n";
	//for (auto &elem : solution)
	//	std::cout << elem << std::endl;
	//return 0;
	
	// fields
	std::vector<double> vel(Nx - 1, 2.);
	std::vector<double> field(Nx, 0.);
	std::vector<double> field_new(Nx, 0.);
	std::vector<double> reconstructed_field(Nx, 0.);
	std::vector<double> time_deriv_GL_1malpha_before(Nx, 0.);
	std::vector<double> time_deriv_GL_1malpha_after(Nx, 0.);
	std::vector<double> sources(Nx, 0.);
	std::vector<double> dphi_dx(Nx);
	for (int ix = 0; ix < Nx; ++ix)
	{
		dphi_dx[ix] = get_initial_field_derivative(dx * ix);
	}
	
	// initial field values: w === 0 initally!
	//get_initial_field_vector(dx, Nx, field);
	//timelayer_fields[0] = field;
	
	// solver
	TVD_scheme solving_scheme(Nx, dx);
	
	// time counters
	double time = 0.;
	int it = 1;
	
	// Grunwald-Letnikov derivative coefficients
	std::vector<double> GL_coeffs_1malpha;
	make_GL_coeff_vec(GL_coeffs_1malpha, 1 - alpha, time_layers_cnt);
	print_field(std::string("GL_coeffs_1malpha.txt"), GL_coeffs_1malpha);
	double dt_1malpha = pow(dt, 1 - alpha);
	double dt_alpha = pow(dt, alpha);
	
	// right-hand part (source and fractional time derivative)
	std::vector<double> rpart(Nx, 0.);
	
	// time cycle
	while(time < finish_time)
	{
		const int N_calc_time_layers = std::min(time_layers_cnt, it + 1);
		
		// fill fractional time derivatives (!!) without first term (!!)
//#pragma omp parallel for
		for(int ix = 0; ix < Nx; ix++)
		{
			time_deriv_GL_1malpha_before[ix] = 0.;
			for (int time_layer = 1; time_layer < N_calc_time_layers; time_layer++)
			{
				time_deriv_GL_1malpha_before[ix] += timelayer_fields[(it - time_layer) % time_layers_cnt][ix] * GL_coeffs_1malpha[time_layer] / dt_1malpha;
			}
		}
		
		// fill sources vector
		get_source_vector(time, dx, Nx, sources);
		
		// fill right part
		// sources
//#pragma omp parallel for
		for (int ix = 0; ix < Nx; ix++)
		{
			// TODO vel[0] --> vel[ix]
			rpart[ix] = sources[ix] - vel[0] * dphi_dx[ix] * pow(time + dt, alpha - 1.) / gamma_alpha;
		}
		// convection part for border
		double time_deriv_GL_border = vel[0] / pow(dt, alpha)
			* (get_border_value(time + dt) - get_initial_field(0.) * pow(time + dt, alpha - 1.) / gamma_alpha);
		
		// solve
		solving_scheme.solve_transfer(vel, time_deriv_GL_1malpha_before, time_deriv_GL_border, field, field_new, dt, alpha, rpart, it);
		
		// save last field
//#pragma omp parallel for
		for (int ix = 0; ix < Nx; ix++)
		{
			timelayer_fields[it % time_layers_cnt][ix] = field_new[ix];
			field[ix] = field_new[ix];
		}
		
		// reconstruction field
		// fill fractional time derivatives
//#pragma omp parallel for
		for (int ix = 0; ix < Nx; ix++)
		{
			time_deriv_GL_1malpha_after[ix] = 0.;
			for (int time_layer = 0; time_layer < N_calc_time_layers; time_layer++)
			{
				time_deriv_GL_1malpha_after[ix] += timelayer_fields[(it - time_layer) % time_layers_cnt][ix] * GL_coeffs_1malpha[time_layer] / dt_1malpha;
			}
		}
		std::vector<double> reconstruction_secondpart(Nx, 0.);
		for (int ix = 0; ix < Nx; ++ix)
		{
			reconstruction_secondpart[ix] = get_initial_field(dx * ix) * pow(time + dt, alpha - 1.) / gamma_alpha;
			reconstructed_field[ix] = time_deriv_GL_1malpha_after[ix] + get_initial_field(dx * ix) * pow(time + dt, alpha - 1.) / gamma_alpha;
		}
		
		// file print
		std::ostringstream stringStream;
		stringStream << "field_" << std::setfill('0') << std::setw(4) << it << ".txt";
		print_field(stringStream.str(), reconstructed_field);
		std::ostringstream stringStream1;
		stringStream1 << "GLfield_" << std::setfill('0') << std::setw(4) << it << ".txt";
		print_field(stringStream1.str(), time_deriv_GL_1malpha_after);
		std::ostringstream stringStream2;
		stringStream2 << "rec2_" << std::setfill('0') << std::setw(4) << it << ".txt";
		print_field(stringStream2.str(), reconstruction_secondpart);

		
		// time counters ++
		time += dt;
		it++;
		std::cout << "Calculated " << time << " secs" << std::endl;
	}
	
	// print final result
	std::string filename("out.txt");
	print_field(filename, reconstructed_field);
	
	return 0;
}

double get_Gamma(const double alpha)
{
	// define alpha1 belonging between 1 and 2
	double alpha1 = alpha - double(int(alpha)) + 1.;
	
	double gamma = 1.;
	if (alpha < 1.)
		gamma /= alpha;
	else if (alpha > 2.)
	{
		double a = alpha;
		while (a > 2.)
		{
			gamma *= a - 1.;
			a -= 1.;
		}
	}
	
	// first term (sum)
	double term1 = 0.;
	// number of terms in first sum
	const int n = 10;
	// threshold for ins=tegral disecting
	const double x0 = 2.;
	for (int i = 0; i < n; ++i)
	{
		double mult = 1. / x0;
		for (int i1 = 0; i1 < i + 1; ++i1)
		{
			mult *= x0 / (alpha1 + i1);
		}
		
		term1 += mult;
	}
	term1 *= pow(x0, alpha1) * exp(-x0);
	
	// second term (integral)
	double term2 = 0.;
	// limit for integral
	const double x_n = 20.;
	// integral step
	const double dx = 0.1;
	double x = x0 + dx * 0.5;
	while (x < x_n)
	{
		term2 += exp(-x) * pow(x, alpha1) / x * dx;
		x += dx;
	}
	
	return gamma * (term1 + term2);
}

void print_field(std::string filename, std::vector<double> &field)
{
	std::ofstream output_file;
	output_file.open(filename, std::ios::out);
	
	for(auto data : field)
	{
		output_file << data << std::endl;
	}
	
	output_file.close();
	
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

void get_initial_field_vector(const double dx, const int Nx, std::vector<double> &field)
{
	for (int ix = 0; ix < Nx; ix++)
	{
		field[ix] = get_initial_field(dx * ix);
	}
	
	return;
}

void get_source_vector(const double time, const double dx, const int Nx, std::vector<double> &source)
{
	for (int ix = 0; ix < Nx; ix++)
	{
		source[ix] = get_source(dx * ix, time);
	}
	
	return;
}