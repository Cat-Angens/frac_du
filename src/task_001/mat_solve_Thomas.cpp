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
	
	bool equal_diag_lens = false;
	assert(n != 0);
	assert(mat.size() == 3);
	assert(mat[0].size() > 2);
	if (mat[0].size() == mat[1].size() - 1 && mat[2].size() == mat[1].size() - 1)
	{
		equal_diag_lens = false;
	}
	else if (mat[0].size() == mat[1].size() && mat[2].size() == mat[1].size())
	{
		equal_diag_lens = true;
		assert(mat[0][0] == 0.);
		assert(mat[2].back() == 0.);
	}
	else
	{
		assert(false);
	}
	
	assert(right_part.size() >= mat[1].size());
	assert(solution.size() >= mat[1].size());
	
	size_t size = mat[1].size();
	std::vector<double> c(size);
	std::vector<double> d(size);
	
	c[0] = mat[2][0] / mat[1][0];
	d[0] = right_part[0] / mat[1][0];
	for (size_t i = 1; i < size; ++i)
	{
		size_t i0, i2;
		if (equal_diag_lens)
		{
			i0 = i;
			i2 = i;
		}
		else
		{
			i0 = i - 1;
			i2 = i;
		}
		
		double denominator = mat[1][i] - mat[0][i0] * c[i - 1];
		if (i != size - 1)
			c[i] = mat[2][i2] / denominator;
		d[i] = (right_part[i] - mat[0][i0] * d[i - 1]) / denominator;
	}
	
	solution.back() = d.back();
	for (int i = (int)size - 2; i >= 0; --i)
	{
		solution[i] = d[i] - c[i] * solution[i + 1];
	}
	
	return;
}