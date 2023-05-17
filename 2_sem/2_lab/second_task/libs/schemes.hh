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
	double T = 0.;
	double X = 0.;
	double t = 0.;
	double x = 0.;
	double *data = nullptr;
	utype u;
    double a = 0.;
    double b = 0.;
    size_t type;

public:
	matrix_t(
        const double T,
        const double X, 
        const int m, 
        const int n, 
        const utype u, 
        const double a, 
        const double b, 
        const int type) :
		    n(n), 
            m(m), 
            T(T), 
            X(X), 
            t((double)T/n), 
            x((double)X/m), 
            data(new double [n * m]), 
            u(u), 
            a(a), 
            b(b),
            type(type) {}

	double &operator()(const size_t i, const size_t j) const noexcept { return data[i*m + j]; }
    double Q(const double _t) const noexcept { return std::pow(_t, a); }
    double K(const double _t) const noexcept { return std::pow(_t, b); }
	void compute();
	void save_img(const std::string, const std::string gif_name);

    double scheme_1(const size_t i, const size_t j, const double c);

	~matrix_t() { delete [] data; }
};

void get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u) noexcept;
};

#endif
