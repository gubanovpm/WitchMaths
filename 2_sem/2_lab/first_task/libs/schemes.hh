#ifndef __schemes_hh__
#define __schemes_hh__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
	int is_div = 0;

	double Lux_scheme_nodiv(const size_t i, const size_t j, const double c) const;
	double Lux_scheme_div(const size_t i, const size_t j, const double c) const;

	double Lux_Wendorf_scheme_nodiv(const size_t i, const size_t j, const double c) const;
	double Lux_Wendorf_scheme_div(const size_t i, const size_t j, const double c) const;
	
	double right_corner_scheme_nodiv(const size_t i, const size_t j, const double c) const;
	double right_corner_scheme_div(const size_t i, const size_t j, const double c) const;
	
	double left_corner_scheme_nodiv(const size_t i, const size_t j, const double c) const;
	double left_corner_scheme_div(const size_t i, const size_t j, const double c) const;
public:
	matrix_t(const double X, const double T,const size_t n, const size_t m, const size_t type, const utype u, const int is_div) :
		n(n), m(m), X(X), T(T), x((double)X/m), t((double)T/n), data(new double [n * m]), u(u), type(type), is_div(is_div) {}

	double &operator()(const size_t i, const size_t j) const noexcept { return data[i*m + j]; }	
	void compute();
	void save_img(const std::string, const std::string gif_name);

	~matrix_t() { delete [] data; }
};

void get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u) noexcept;
};

#endif
