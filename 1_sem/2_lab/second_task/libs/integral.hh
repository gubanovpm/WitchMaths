#ifndef __integral_hh__
#define __integral_hh__

#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

namespace witch_math {

using function_type = std::function<double(const double)>;

struct IIntegral {
public:
	IIntegral (const function_type &f) : _f(f) {}	
	virtual ~IIntegral() {}

	double value (const std::vector<double> &x) const {
		double result = 0;	
		for (size_t k = 0; k < x.size() - 1; ++k)
			result += point_value(x[k], x[k+1] - x[k]);
		return result;
	}

	void set_function(const function_type f) { _f = f; }
protected:
	function_type _f;
	virtual double point_value (const double x, const double h) const = 0;
};

struct rectangle_integral final : public IIntegral {
public:
	rectangle_integral(const function_type &f) : IIntegral(f) {}
private:
	double point_value (const double x, const double h) const override {	
		double res = _f(x + h/2);
		return h * res;
	}
};

struct trapezoid_integral final : public IIntegral {
public:
	trapezoid_integral(const function_type &f) : IIntegral(f) {}
private:
	double point_value (const double x, const double h) const override {
		return h/2 * (_f(x) + _f(x+h));
	}
};

struct simpsons_integral final : public IIntegral {
public:
	simpsons_integral(const function_type &f) : IIntegral(f) {}
private:
	double point_value (const double x, const double h) const override {
		return h/6 * (_f(x) + 4*_f(x+h/2) + _f(x+h));
	}
};

struct three_eighths_integral final : public IIntegral {
public:
	three_eighths_integral(const function_type &f) : IIntegral(f) {}
private:
	double point_value (const double x, const double h) const override {
		return h/8*(_f(x) + 3*_f(x+h/3) + 3*_f(x+2*h/3) + _f(x+h));
	}
};

}

#endif
