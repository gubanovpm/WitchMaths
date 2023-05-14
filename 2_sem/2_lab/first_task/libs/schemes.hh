#ifndef __schemes_hh__
#define __schemes_hh__

#include <iostream>
#include <matplot/matplot.h>
#include <pthread.h>

#ifndef THREAD_COUNT
#define THREAD_COUNT 6
#endif

using utype = std::function<double(const double)>;

namespace WitchMath {

enum SchemeT {
	Lux = 0, 
	LuxWendorf, 
	RightCorner, 
	LeftCorner
};

struct matrix_t final{
private:
	size_t __n = 0;
	size_t __m = 0;
	double *__data = nullptr;
	utype __u;

	double Lux_scheme(const size_t i, const size_t j) const noexcept;
	double Lux_Wendorf_scheme(const size_t i, const size_t j) const noexcept;
	double right_corner_scheme(const size_t i, const size_t j) const noexcept;
	double left_corner_scheme(const size_t i, const size_t j) const noexcept;
public:
	matrix_t(const size_t n, const size_t m, const utype u) :
		__n(n), __m(m), __data(new double [n * m]), __u(u) {}

	double &operator()(const size_t i, const size_t j) const noexcept { return __data[i*__m + j]; }	
	size_t getn() const noexcept { return __n; }
	size_t getm() const noexcept { return __m; }
	double u(const double x) const {return __u(x);}

	~matrix_t() { delete [] __data; }
};

struct thread_arg_t final {
	size_t id;
	size_t np;
	size_t tp;
	matrix_t *matrix;
};

void *compute(void *thread_args);

};

#endif
