#ifndef __seidel_hh__
#define __seidel_hh__

#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <matplot/matplot.h>

using namespace boost::numeric::ublas;

matrix<double> get_right_matrix(size_t n, double a);
matrix<double> get_right_vector(size_t n, double a);
matrix<double> get_right_answer(size_t n, double a);

struct norm {
    double euclid(const matrix<double> &) const;
};

struct seidel_method {
private:
    matrix<double> _c;
    matrix<double> _d;
    std::vector<matrix<double>> _itter = {};

    void get_new_k_itter (const size_t k);
public:
    seidel_method(const matrix<double> &A, const matrix<double> &f);
    matrix<double> &get_k_itter (const size_t k);
};

#endif