///////////////////////////
// Solver of Linear algebraic 3-diagonal system by Thomas method

class mat_solve_Thomas;

#pragma once

#include <vector>

class mat_solve_Thomas
{
public:
	mat_solve_Thomas(size_t n_ = 0);
	~mat_solve_Thomas();
private:
	size_t n; // matrix dimension
	
public:
	
	void set_dimension(size_t n_);
	// Решение СЛАУ
	void solve(
		const std::vector<double> &mat,
		const std::vector<double> &right_part,
		std::vector<double> &solution
		);
	
	void solve(
		const std::vector<std::vector<double>> &mat,
		const std::vector<double> &right_part,
		std::vector<double> &solution
		);
};
