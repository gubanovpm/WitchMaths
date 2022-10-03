#include "../libs/seidel.hh"

matrix<double> get_right_matrix(size_t n, double a) {
    matrix<double> A(n, n);
    for (size_t i = 0; i < A.size1(); ++i)
        for (size_t j = 0; j < A.size2(); ++j) {
            if      (i == j)      A(i, j) = 2;
            else if (j == i + 1)  A(i, j) = -1 - a;
            else if (j == i - 1)  A(i, j) = -1 + a;
            else                  A(i, j) = 0;
        }
    return A;
}

matrix<double> get_right_vector(size_t n, double a) {
    matrix<double> f(n, 1);
    if   (n == 1) 
        f(0, 0) = 1 - a;
    else {
        f(0, 0)     = 1 - a;
        f(n - 1, 0) = 1 + a;
    }
    for (size_t i = 1; i + 1 < f.size1(); ++i)
        f(i, 0) = 0;
    return f;
}

matrix<double> get_right_answer(size_t n, double a) {
    matrix<double> answer(n, 1);
    if (n == 1) { answer(0, 0) = (1 - a) / 2; return answer; }
    for (size_t i = 0; i < n; ++i) answer(i, 0) = 1;
    return answer;
}

seidel_method::seidel_method(const matrix<double> &A, const matrix<double> &f) {
    _c = matrix<double> {A.size1(), A.size2()};
    _d = matrix<double> {f.size1(), f.size2()};

    for (size_t i = 0; i < _c.size1(); ++i) {
        for (size_t j = 0; j < _c.size2(); ++j) {
            _c(i, j) = ((i != j) ? -A(i, j)/A(i, i) : 0 );
        }
        _d(i, 0) = f(i, 0) / A(i, i);
    }
    
    _itter.push_back(_d);
}

matrix<double> &seidel_method::get_k_itter(const size_t k) {
    if (k < _itter.size()) return _itter[k];
    
    get_new_k_itter(k);    
    return _itter[k];
}

void seidel_method::get_new_k_itter (const size_t k) {
    matrix<double> new_itter = _d;
    for (size_t i = 0; i < new_itter.size1(); ++i) {
        for (size_t j = 0; (i != 0) && j < i - 1; ++j) 
            new_itter(i, 0) += _c(i, j) * new_itter(j, 0);
        for (size_t j = 0; j < new_itter.size1(); ++j) 
            new_itter(i, 0) += _c(i, j) * get_k_itter(k - 1)(j, 0);
    }
    _itter.push_back(new_itter);
}

double norm::euclid(const matrix<double> &A) const {
    double result = 0;
    if (A.size2() == 1) 
        for (size_t i = 0; i < A.size1(); ++i)
            result += A(i, 0) * A(i, 0);
    return std::sqrt(result);
}