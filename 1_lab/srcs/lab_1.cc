#include "../libs/lab_1.hh"

size_t factorial(const  size_t n) {
    size_t res = 1;
    for (size_t i = 2; i <= n; ++i) res *= i;
    return res;
}
size_t C(const size_t n,const size_t k) {
    return factorial(n) / factorial(k) / factorial(n - k);
}

double central_derivative(const std::function<double(const double)> &f, const size_t n, const double x, const double dt) {
    double cur_sum = 0;
    if (n == 0) return f(x);
    for (size_t i = 0; i <= n; ++i) {
        cur_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * f(x + (n - double(2*i)) * dt);
    }
    return cur_sum / std::pow(2*dt, n);
}

double forward_derivative(const std::function<double(const double)> &f, const size_t n, const double x, const double dt) {
    double cur_sum = 0;
    if (n == 0) return f(x);
    for (size_t i = 0; i <= n; ++i) {
        cur_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * f(x + (n - double(i)) * dt);
    }
    return cur_sum / std::pow(dt, n);
}

double back_derivative(const std::function<double(const double)> &f, const size_t n, const double x, const double dt) {
    double cur_sum = 0;
    if (n == 0) return f(x);
    for (size_t i = 0; i <= n; ++i) {
        cur_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * f(x - (double(i)) * dt);
    }
    return cur_sum / std::pow(dt, n);
}

double maclaurin_series(const std::function<double(const double)> &f, 
                        const std::function<double(const std::function<double(const double)> &, const size_t, const double, const double)> &der_f,
                        double x, size_t n, const double dt) {
    double cur_sum = 0;
    for (size_t i = 0; i <= n; ++i) {
        cur_sum += double(std::pow(x, i)) * der_f(f, i, 0, dt) / factorial(i);
    }
    return cur_sum;
}

size_t get_better_approx_n(const std::function<double(const double)> &f, 
                           const std::function<double(const std::function<double(const double)> &, const size_t, const double, const double)> &der_f,
                           const size_t kmax, const float dt, const float begin, const float end) {
    bool is_in_epsilon = 0;
    size_t n = 0;
    for (size_t k = 1; k <= kmax; ++k) {
        is_in_epsilon = 1;
        for (float cur = begin; cur <= end; cur += dt) {
            if (std::abs(maclaurin_series(f, der_f, cur, k, dt) - f(cur)) > 2 * dt) {
                is_in_epsilon = 0; break;
            }
        }
        if (is_in_epsilon) n = k;
    }
    return n;
}