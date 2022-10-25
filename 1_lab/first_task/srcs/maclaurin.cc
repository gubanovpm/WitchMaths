#include "../libs/maclaurin.hh"

size_t factorial(const  size_t n) {
    size_t res = 1;
    for (size_t i = 2; i <= n; ++i) res *= i;
    return res;
}
size_t C(const size_t n,const size_t k) {
    return factorial(n) / factorial(k) / factorial(n - k);
}

size_t IDerivative::__find_it_key(const double dh){
    for (size_t i = 0; i < _val.size(); ++i)
        if (std::abs(_val[i].first - dh) <= EPSILON) { return i; }
    return -1;
}

double c_derivative::value(const size_t n, const arg_type x, const double dh) {
    size_t key = __find_it_key(dh);
    if (key != -1 && _val[key].second.size() > n) return _val[key].second[n];
    if (key == -1) { key = _val.size() ; _val.push_back({dh, {}}); }

    if (_val[key].second.size() == 0) {
        _val[key].second[0] = _f(x); return _val[key].second[0];
    }

    double t_sum = 0;
    for (size_t i = 0; i <= n; ++i) {
        t_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * _f(x + (n - double(2*i)) * dh);
    }

    _val[key].second[n] = t_sum / std::pow(2*dh, n);
    return _val[key].second[n];
}

double b_derivative::value(const size_t n, const arg_type x, const double dh) {
    size_t key = __find_it_key(dh);
    if (key != -1 && _val[key].second.size() > n) return _val[key].second[n];
    if (key == -1) { 
        key = _val.size(); _val.push_back({dh, {}});
    }

    if (_val[key].second.size() == 0) {
        _val[key].second[0] = _f(x);
        return _val[key].second[0];
    }

    double t_sum = 0;
    for (size_t i = 0; i <= n; ++i) {
        t_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * _f(x - (double(i)) * dh);
    }
    _val[key].second[n] = t_sum / std::pow(dh, n);
    return _val[key].second[n];
}

double f_derivative::value(const size_t n, const arg_type x, const double dh) {
    size_t key = __find_it_key(dh);
    if (key != -1 && _val[key].second.size() > n) return _val[key].second[n];
    if (key == -1) { 
        key = _val.size(); _val.push_back({dh, {}});
    }

    if (_val[key].second.size() == 0) {
        _val[key].second[0] = _f(x);
        return _val[key].second[0];
    }

    double t_sum = 0;
    for (size_t i = 0; i <= n; ++i) {
        t_sum += double(C(n, i)) * ((i % 2) ? -1 : 1) * _f(x + (n - double(i)) * dh);
    }
    _val[key].second[n] = t_sum / std::pow(dh, n);
    return _val[key].second[n];
}

maclaurin_series::maclaurin_series(const function_type &f, const D_TYPE type) {
    switch (type) {
        case CENTRAL: { _der = new c_derivative(f); break; }
        case FORWARD: { _der = new f_derivative(f); break; }
        case BACK   : { _der = new b_derivative(f); break; }
        default: throw;
    }
}

size_t maclaurin_series::__find_key(const double x) {
    for (size_t i = 0; i < _val.size(); ++i)
        if (std::abs(_val[i].first - x) <= EPSILON) { return i; }
    return -1;
}

double maclaurin_series::value(const size_t n, const double x, const double dt) {
    size_t key = __find_key(x);
    if (key != -1 && _val[key].second.size() > n) return _val[key].second[n];
    if (key == -1) { key = _val.size(); _val.push_back({x, {}}); }

    if (_val[key].second.size() == 0) _val[key].second.push_back(_der->value(0, 0, dt));
    for (size_t k = _val[key].second.size(); k <= n; ++k) {
        _val[key].second.push_back(_val[key].second[k - 1] + std::pow(x, k) * _der->value(k, 0, dt) / factorial(k));
    }
    return _val[key].second[n];

}