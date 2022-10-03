#ifndef __maclaurin_hh__
#define __maclaurin_hh__

#include <iostream>
#include <functional>
#include <cmath>

#include <matplot/matplot.h>

double central_derivative(const std::function<double(const double)> &, const size_t n, const double x, const double dt);
double forward_derivative(const std::function<double(const double)> &, const size_t n, const double x, const double dt);
double back_derivative(const std::function<double(const double)> &, const size_t n, const double x, const double dt);

double maclaurin_series(const std::function<double(const double)> &f, 
                        const std::function<double(const std::function<double(const double)> &, const size_t, const double, const double)> &der_f,
                        double x, size_t n, const double dt);

size_t get_better_approx_n(const std::function<double(const double)> &f, 
                           const std::function<double(const std::function<double(const double)> &, const size_t, const double, const double)> &der_f,
                           const size_t kmax, const float dt, const float begin, const float end);

size_t factorial(const size_t);
size_t C(const size_t,const size_t);

#endif