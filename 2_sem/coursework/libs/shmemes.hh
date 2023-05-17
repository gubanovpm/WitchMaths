#ifndef __schmemes_hh__
#define __schmemes_hh__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <matplot/matplot.h>
#include <utility>

const double R0  = 100;
const double T0  = 0.3;
const double E0  = 3.15e9;
const double ro0 = 5.45e-4;
const double P0  = 1.15e6;
const double u0  = 0;
const double g   = 1.4;

using utype = std::function<double(const double)>;

namespace WitchMath {

struct matrix_t final{
private:
	size_t n = 0;
	size_t m = 0;
	double T = 0.;
	double X = 0.;
	double tau = 0.;
	double h = 0.;
	std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> E;
    std::vector<std::vector<double>> P;
    std::vector<std::vector<double>> ro;

public:
	matrix_t(
        const double T,
        const double X, 
        const int n, 
        const int m
        ) :
		    n(n),
            m(m),
            T(T),
            X(X),
            tau((double)T/n),
            h((double)X/m),
            u(n),
            E(n),
            P(n),
            ro(n) {
                for (size_t i = 0; i < n; ++i) { 
                    u [i].resize(m+1),
                    E [i].resize(m+1),
                    P [i].resize(m+1),
                    ro[i].resize(m+1);
                }
                for (size_t j = 0; j <= m; ++j) {
                    u [0][j] = u0,
                    E [0][j] = ((j*h <= R0) ? E0  : 0),
                    P [0][j] = ((j*h <= R0) ? P0  : 0),
                    ro[0][j] = ((j*h <= R0) ? ro0 : 0);
                }
            }

	void compute();
	void view_plot(const size_t num) noexcept;

    void scheme_u (const size_t i, const size_t j, const double c) noexcept;
    void scheme_u (const size_t i, const double c) noexcept;

    void scheme_ro(const size_t i, const size_t j, const double c) noexcept;
    void scheme_ro(const size_t i, const double c) noexcept;

    void scheme_EP(const size_t i, const size_t j, const double c) noexcept;
    void scheme_EP(const size_t i, const double c) noexcept;
	~matrix_t() {}
};

void get_args(const int argc, const char *argv[], matrix_t **matrix) noexcept;
void check_curant(const double c, const double tau, const double h) noexcept;
};

#endif
