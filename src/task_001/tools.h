#include <string>
#include <vector>

void fprint_vector(std::string filename, const std::vector<double> &field);

void fprint_vector(const std::string prefix, int iter, const std::vector<double> &field);