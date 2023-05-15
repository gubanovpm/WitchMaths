#ifndef __schemes_hh__
#define __schemes_hh__

#include <iostream>
#include <matplot/matplot.h>

using utype = std::function<double(const double)>;

namespace WitchMath {

struct matrix_t final{
private:
	size_t n = 0;
	size_t m = 0;
	double X = 0.;
	double T = 0.;
	double x = 0.;
	double t = 0.;
	double *data = nullptr;
	utype u;
	size_t type = 0;

	double Lux_scheme(const size_t i, const size_t j, const double c) const;
	double Lux_Wendorf_scheme(const size_t i, const size_t j) const;
	double right_corner_scheme(const size_t i, const size_t j) const;
	double left_corner_scheme(const size_t i, const size_t j) const;
public:
	matrix_t(const double X, const double T,const size_t n, const size_t m, const size_t type, const utype u) :
		n(n), m(m), X(X), T(T), x((double)X/m), t((double)T/n), data(new double [n * m]), u(u), type(type) {}

	double &operator()(const size_t i, const size_t j) const noexcept { return data[i*m + j]; }	
	void compute();
	void save_img(const std::string);

	~matrix_t() { delete [] data; }
};

void get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u) noexcept;
};

#endif
