
#include <vector>

class TVD_scheme;


#pragma once

#include "mat_solve_Thomas.h"

class TVD_scheme
{
	
	const int Nx;
	const double dx;
	
	std::vector<double>  Phi;
	std::vector<std::vector<double>>  mat_F;        // 3-диагональная матрица конвекции
	std::vector<std::vector<double>>  mat_full;     // полная 3-диагональная матрица системы
	std::vector<double>               right_part;   // вектор правой части
	
	mat_solve_Thomas                  *solver;      // решатель СЛАУ методом прогонки
	
	double get_r_i_biased(const std::vector<double> &field, const std::vector<double> &vel_edg, int ix) const;
	void fill_F(const std::vector<double> &field, const std::vector<double> &vel_edg);
	void fill_F_without_tvd(const std::vector<double> &field, const std::vector<double> &vel_edg);
	void fill_rpart(const std::vector<double> &field_GL_derivative, const std::vector<double>& sources, const double GL_derivative_border);
	void fill_full_matrix(const double dt_1malpha);
	void fill_Phi(const std::vector<double> &field, const std::vector<double> &vel_edg);
	double tvd_limit(const double r) const;
	
public:
	
	TVD_scheme(int nx, double dx);
	
	void solve_transfer_explicitly(const std::vector<double> &vel, const std::vector<double> &field_GL_derivative, const double GL_derivative_border, const std::vector<double> & field_old, std::vector<double> & field_new, const double dt, const double alpha, const std::vector<double> & right_part);
	
};



