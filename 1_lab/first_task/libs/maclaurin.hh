#ifndef __maclaurin_hh__
#define __maclaurin_hh__

#include <iostream>
#include <functional>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <matplot/matplot.h>

double central_derivative(const std::function<double(const double)> &, const size_t n, const double x, const double dt);
double forward_derivative(const std::function<double(const double)> &, const size_t n, const double x, const double dt);
double back_derivative   (const std::function<double(const double)> &, const size_t n, const double x, const double dt);

// for the future
// using arg_type      = double;
// using function_type = std::function<double(const arg_type)>;

// struct derivative {
// public:
//     enum DERIVATIVE_T {CENTRAL, BACK, FORWARD, D_TYPES};

//     std::unordered_map<double, std::vector<double>> _val[D_TYPES];

//     derivative(const function_type &f): _f(f) {}
//     double value(const size_t n, const arg_type x, const double dh, const DERIVATIVE_T type); 
// private:
//     const function_type &_f;

//     double __c(const size_t n, const arg_type x, const double dh);
//     double __b(const size_t n, const arg_type x, const double dh);
//     double __f(const size_t n, const arg_type x, const double dh);
// };

// struct maclaurin_series final {
// private:
//     const function_type &_f;
//     std::vector<double> _val[derivative::D_TYPES] = {};
// public:
//     maclaurin_series(const function_type &f) : _f(f) {}

//     double get_n(const size_t n, const double x, const double dt);
//     double value(const size_t n, const double x, const double dt); 
// };
double maclaurin_series(const std::function<double(const double)> &f, 
                        const std::function<double(const std::function<double(const double)> &, const size_t, const double, const double)> &der_f,
                        double x, size_t n, const double dt);

size_t get_better_approx_n(const std::function<double(const double)> &f, 
                           const std::function<double(const std::function<double(const double)> &, const size_t, const double, const double)> &der_f,
                           const size_t kmax, const float dt, const float begin, const float end);

size_t factorial(const size_t);
size_t C(const size_t,const size_t);

#endif