
#include <vector>

class TVD_scheme;


#pragma once

class TVD_scheme
{
	
	const int Nx;
	const double dx;
	
	std::vector<double>  Phi;
	std::vector<std::vector<double>>  mat_F;     // 3-диагональная матрица
	
	double get_r_i(const std::vector<double> &field, int ix, bool positive_stream) const;
	double get_r_i_biased(const std::vector<double> &field, const std::vector<double> &vel_edg, int ix) const;
	void fill_F(const std::vector<double> &field, const std::vector<double> &vel_edg);
	void fill_Phi(const std::vector<double> &field, const std::vector<double> &vel_edg);
	double tvd_limit(const double r) const;
	
public:
	
	TVD_scheme(int nx, double dx);
	
	void solve_transfer_explicitly(const std::vector<double> &vel, const std::vector<double> &field_old, std::vector<double> &field_new, const double dt);
	
};



