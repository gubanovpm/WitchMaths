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
	std::vector<std::vector<double>> data;
	utype u;
    utype g;
    double a = 0.;
    double b = 0.;

public:
	matrix_t(
        const double T,
        const double X, 
        const int m, 
        const int n, 
        const utype u,
        const utype g,
        const double a, 
        const double b) :
		    n(n), 
            m(m), 
            T(T), 
            X(X), 
            t((double)T/n), 
            x((double)X/m), 
            data(n), 
            u(u), 
            g(g),
            a(a), 
            b(b) {
                for (auto &line: data) line.resize(m); 
            }

    double Q(const double _t) const noexcept { return std::pow(_t, a); }
    double K(const double _t) const noexcept { return std::pow(_t, b); }
    double amn(const size_t i, const size_t j) const noexcept { 
        size_t jmx = ((j+1 == m) ? j : j + 1);
        return 2 * K(data[i][j]) * K(data[i][jmx]) / (K(data[i][j]) + K(data[i][jmx])); 
    }
	void compute();
	void save_img() noexcept;

    double scheme_1(const size_t i, const size_t j, const double c);

	~matrix_t() { }
};

void get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u, const utype g) noexcept;
};

#endif
