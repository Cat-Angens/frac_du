
#include "TVD_scheme.h"
#include <fstream>

TVD_scheme::TVD_scheme(int nx, double dx)
	: Nx(nx),
	dx(dx)
{
	mat_F.resize(Nx, std::vector<double>(3, 0.));
	Phi.resize(Nx - 1, 0.);
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
			mat_F[ix    ][1] += - phi_term;
			mat_F[ix    ][2] +=   phi_term + vel_xy;
			
			mat_F[ix + 1][0] +=   phi_term;
			mat_F[ix + 1][1] += - phi_term - vel_xy;
		}
		else
		{
			mat_F[ix    ][1] += - phi_term + vel_xy;
			mat_F[ix    ][2] +=   phi_term;
			
			mat_F[ix + 1][0] +=   phi_term - vel_xy;
			mat_F[ix + 1][1] += - phi_term;
		}
	}
	
	return;
}

double TVD_scheme::tvd_limit(const double r) const
{
	if (r < 0)
		return 0.;
	
	return fmax(fmin(1.0, 2 * r), fmin(2.0, r));
	
}

void TVD_scheme::solve_transfer_explicitly(const std::vector<double> &vel, const std::vector<double> &field_old, std::vector<double> &field_new, const double dt, const std::vector<double> &right_part)
{
	fill_F(field_old, vel);
	
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
#endif
	
#pragma omp parallel for
	for (int ix = 0; ix < Nx; ix++)
	{
		field_new[ix] = right_part[ix] * dt;
		
		for (int adj = 0; adj < 3; adj++)
		{
			if ((ix == 0 && adj == 0) || (ix == Nx - 1 && adj == 2))
			{
				continue;
			}
			field_new[ix] -= dt * mat_F[ix][adj] * field_old[ix + adj - 1];
		}
	}
}

void TVD_scheme::fill_Phi(const std::vector<double> &field, const std::vector<double> &vel_edg)
{
#pragma omp parallel for
	for (unsigned int ix = 0; ix < Nx - 1; ix++)
	{
		Phi[ix] = tvd_limit(get_r_i_biased(field, vel_edg, ix));
	}
	
}
