#include "mat_solve_Thomas.h"
#include <assert.h>

mat_solve_Thomas::mat_solve_Thomas(size_t n_) : n(n_)
{

}

mat_solve_Thomas::~mat_solve_Thomas()
{
}

void mat_solve_Thomas::set_dimension(size_t n_)
{
	n = n_;
}

void mat_solve_Thomas::solve(
	const std::vector<double> &mat,
	const std::vector<double> &right_part,
	std::vector<double> &solution
	)
{
	if (n == 0)
		n = (mat.size() + 2) / 3;
	
	assert(n != 0);
	assert(mat.size() % 3 == 1);
	assert(mat.size() > 4);
	
	
	
	
	return;
}

void mat_solve_Thomas::solve(
	const std::vector<std::vector<double>> &mat,
	const std::vector<double> &right_part,
	std::vector<double> &solution
	)
{
	
	if (n == 0)
		n = mat[1].size();
	
	assert(n != 0);
	assert(mat.size() == 3);
	assert(mat[0].size() > 2);
	assert(mat[1].size() == mat[0].size() + 1);
	assert(mat[2].size() == mat[0].size());
	
	
	
	
	return;
}