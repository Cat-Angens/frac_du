
#include "TVD_scheme.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#define printmat

TVD_scheme::TVD_scheme(int nx, double dx)
	: Nx(nx),
	dx(dx)
{
	Phi.resize(Nx - 1, 0.);
	mat_F.resize(Nx, std::vector<double>(3, 0.));
	mat_full.resize(3, std::vector<double>(Nx, 0.));
	right_part.resize(Nx, 0.);
	
	solver = new mat_solve_Thomas(Nx);
	
	std::ofstream output_file;
	output_file.open(std::string("mat_f.txt"), std::ios::out);
	output_file.close();
}

double TVD_scheme::get_r_i_biased(const std::vector<double> &field, const std::vector<double> &vel_edg, int ix) const
{
	double ri_curr;
	double field_up, field_dn, field_cr;
	
	double vel = vel_edg[ix];
	
	if((vel > 0 && ix != 0) || ix == Nx - 2)
	{
		field_cr = field[ix    ];
		field_up = field[ix + 1];
		field_dn = field[ix - 1];
	} else
	{
		field_cr = field[ix + 1];
		field_up = field[ix    ];
		field_dn = field[ix + 2];
	}
	
	// TODO: сравнивать до заданной точности
	if(field_cr == field_up){
		if(field_cr - field_dn < 0.)
			ri_curr = -1e6;
		else if(field_cr - field_dn > 0.)
			ri_curr = 1e6;
		else
			ri_curr = 1.;
	}
	else{
		ri_curr = (field_cr - field_dn) / (field_up - field_cr);
		if(std::isinf(ri_curr) || std::isnan(ri_curr))
			ri_curr = 1;
	}
	return ri_curr;
}

void TVD_scheme::fill_F(const std::vector<double> &field, const std::vector<double> &vel_edg)
{
#pragma omp parallel for
	for(int ix = 0; ix < Nx; ix++)
	{
		for (int adj = 0; adj < 3; adj++)
		{
			mat_F[ix][adj] = 0;
		}
	}
	
	fill_Phi(field, vel_edg);
	
#pragma omp parallel for
	// расчет переноса через ребра выбранного направления
	for(int ix = 0; ix < Nx - 1; ix++)
	{
		// рассматривается взаимодействие между ячейками ix и ix + 1
		
		double vel_xy = vel_edg[ix] / dx;
		
		double phi_term = 0.5 * fabs(vel_xy) * Phi[ix];
		
		if(vel_xy < 0)
		{
			mat_F[ix    ][1] -= - phi_term;
			mat_F[ix    ][2] -=   phi_term + vel_xy;
			
			mat_F[ix + 1][0] -=   phi_term;
			mat_F[ix + 1][1] -= - phi_term - vel_xy;
		}
		else
		{
			mat_F[ix    ][1] -= - phi_term + vel_xy;
			mat_F[ix    ][2] -=   phi_term;
			
			mat_F[ix + 1][0] -=   phi_term - vel_xy;
			mat_F[ix + 1][1] -= - phi_term;
		}
	}
	
	return;
}

void TVD_scheme::fill_F_without_tvd(const std::vector<double> &field, const std::vector<double> &vel_edg)
{
#pragma omp parallel for
	for(int ix = 0; ix < Nx; ix++)
	{
		for (int adj = 0; adj < 3; adj++)
		{
			mat_F[ix][adj] = 0;
		}
	}
	
	// ix == 0      <-> ix == 1
	mat_F[0][2] += vel_edg[0] / dx;
	mat_F[1][0] -= 0.5 * vel_edg[0] / dx;
	
	// расчет переноса через ребра выбранного направления
#pragma omp parallel for
	for(int ix = 1; ix < Nx - 1; ix++)
	{
		// рассматривается взаимодействие между ячейками ix и ix + 1
		
		double vel_xy = 0.5 * vel_edg[ix] / dx;
		
		mat_F[ix    ][2] += vel_xy;
		mat_F[ix + 1][0] -= vel_xy;
	}
	
	// ix == Nx - 2 <-> ix == Nx - 1
	mat_F[Nx - 1][0] = -vel_edg[Nx - 2] / dx;
	mat_F[Nx - 1][1] =  vel_edg[Nx - 2] / dx;
	
	return;
}

void TVD_scheme::fill_rpart(const std::vector<double> & field_old, const std::vector<double>& field_GL_derivative, const std::vector<double>& sources, const double GL_derivative_border, const double dt)
{
	for (int ix = 1; ix < Nx; ix++)
	{
		right_part[ix] = field_old[ix] + dt * sources[ix];
		for (int adj = 0; adj < 3; ++adj)
		{
			if (ix == Nx - 1 && adj == 2)
				continue;
			right_part[ix] -= dt * mat_F[ix][adj] * field_GL_derivative[ix - 1 + adj];
		}
	}
	right_part[0] = dt * (sources[0] + GL_derivative_border);
}

void TVD_scheme::fill_full_matrix(const double dt_1malpha)
{
	for (int ix = 0; ix < Nx; ix++)
	{
		for (int adj = 0; adj < 3; adj++)
		{
			mat_full[adj][ix] = dt_1malpha * mat_F[ix][adj];
			if (adj == 1)
				mat_full[adj][ix] += 1.;
		}
	}
}

double TVD_scheme::tvd_limit(const double r) const
{
	if (r < 0)
		return 0.;
	
	return fmax(fmin(1.0, 2 * r), fmin(2.0, r));
	
}

void TVD_scheme::solve_transfer_explicitly(
	const std::vector<double> &vel,
	const std::vector<double> &field_GL_derivative,
	const double GL_derivative_border,
	const std::vector<double> & field_old,
	std::vector<double> & field_new,
	const double dt,
	const double alpha,
	const std::vector<double> & sources,
	const int it)
{
	double dt_1malpha = pow(dt, 1 - alpha);
	//fill_F(field_old, vel);
	fill_F_without_tvd(field_old, vel);
	fill_rpart(field_old, field_GL_derivative, sources, GL_derivative_border, dt);
	fill_full_matrix(dt_1malpha);
	
#ifdef printmat
	std::ofstream output_file;
	output_file.open(std::string("mat_f.txt"), std::ios::app);
	
	for(auto row : mat_F)
	{
		for(auto elem : row)
		{
			output_file << elem << "\t";
		}
		output_file << std::endl;
	}
	output_file << "################ FINISHED F ################" << std::endl << std::endl;
	output_file.close();
#endif
	
	solver->solve(mat_full, right_part, field_new);
//#pragma omp parallel for
//	for (int ix = 0; ix < Nx; ix++)
//	{
//		field_new[ix] = field_old[ix] + sources[ix] * dt;
//		
//		for (int adj = 0; adj < 3; adj++)
//		{
//			if ((ix == 0 && adj == 0) || (ix == Nx - 1 && adj == 2))
//			{
//				continue;
//			}
//			field_new[ix] -= dt * mat_F[ix][adj] * field_GL_derivative[ix + adj - 1];
//		}
//	}
	
	std::ostringstream stringStream1;
	stringStream1 << "GL_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream1.str(), field_GL_derivative);
	std::ostringstream stringStream2;
	stringStream2 << "rpart_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream2.str(), right_part);
	std::ostringstream stringStream3;
	stringStream3 << "sources_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream3.str(), sources);
	std::ostringstream stringStream4;
	stringStream4 << "mat0_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream4.str(), mat_full[0]);
	std::ostringstream stringStream5;
	stringStream5 << "mat1_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream5.str(), mat_full[1]);
	std::ostringstream stringStream6;
	stringStream6 << "mat2_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream6.str(), mat_full[2]);
	std::ostringstream stringStream7;
	stringStream7 << "w_" << std::setfill('0') << std::setw(4) << it << ".txt";
	print_field(stringStream7.str(), field_new);
	
}

void TVD_scheme::fill_Phi(const std::vector<double> &field, const std::vector<double> &vel_edg)
{
#pragma omp parallel for
	for (unsigned int ix = 0; ix < Nx - 1; ix++)
	{
		Phi[ix] = tvd_limit(get_r_i_biased(field, vel_edg, ix));
	}
	
}

void TVD_scheme::print_field(std::string filename, const std::vector<double> &field) const
{
	std::ofstream output_file;
	output_file.open(filename, std::ios::out);

	for (auto data : field)
	{
		output_file << data << std::endl;
	}

	return;
}