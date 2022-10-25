#ifndef __maclaurin_hh__
#define __maclaurin_hh__

#include <matplot/matplot.h>
#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#define EPSILON 1e-4

using arg_type      = double;
using function_type = std::function<double(const arg_type)>;

struct IDerivative {
public:
    IDerivative(const function_type &f):_f(f) {};

    virtual double value(const size_t n, const arg_type x, const double dh) = 0; 
    size_t size() { return _val.size(); }
    
    virtual ~IDerivative() {}
protected:
    const function_type &_f;
    std::vector<std::pair<double, std::map<size_t, double>>> _val;
    size_t __find_it_key(const double dh);
};

struct c_derivative final: public IDerivative {
public:
    c_derivative(const function_type &f): IDerivative(f) {}
    double value(const size_t n, const arg_type x, const double dh) override;
};

struct b_derivative final: public IDerivative {
public:
    b_derivative(const function_type &f): IDerivative(f) {}
    double value(const size_t n, const arg_type x, const double dh) override;
};

struct f_derivative final: public IDerivative {
public:
    f_derivative(const function_type &f): IDerivative(f) {}
    double value(const size_t n, const arg_type x, const double dh) override;
};

struct maclaurin_series final {
private:
    IDerivative *_der = nullptr;
    std::vector<std::pair<double, std::vector<double>>> _val = {};

    size_t __find_key(const double x);
public:
    enum D_TYPE {CENTRAL, BACK, FORWARD, D_TYPES};
    maclaurin_series(const function_type &f, const D_TYPE type);
    size_t size() { return _der->size(); }

    double value(const size_t n, const double x, const double dt); 
};

size_t factorial(const size_t);
size_t C(const size_t,const size_t);

#endif