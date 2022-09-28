#include "../libs/lab_1.hh"

using namespace matplot;

#define EXP

int main() {
    #ifdef SIN
        auto f = [](double x){ return std::sin(x); };
    #else 
        auto f = [](double x){ return std::exp(x); };
    #endif
    double begin = 0, end = 1, dt = 0.001;
    size_t kmax = 8;
    
    std::cout << "Better n for central derivative = " << 
                 get_better_approx_n(f, central_derivative, kmax, dt, begin, end) << std::endl;
    std::cout << "Better n for forward derivative = " << 
                 get_better_approx_n(f, forward_derivative, kmax, dt, begin, end) << std::endl;
    std::cout << "Better n for back    derivative = " << 
                 get_better_approx_n(f, back_derivative   , kmax, dt, begin, end) << std::endl;
    
    std::vector<size_t> x   = {};
    std::vector<double> y_1 = {};
    std::vector<double> y_2 = {};
    std::vector<double> y_3 = {};
    for (size_t k = 1; k <= kmax; ++k) {
        double sum_1 = 0, sum_2 = 0, sum_3 = 0; size_t count = 0;
        for (float cur = begin; cur <= end; cur += dt) {
            sum_1 += std::abs(maclaurin_series(f, central_derivative, cur, k, dt) - f(cur));
            if (k != kmax) {
                sum_2 += std::abs(maclaurin_series(f, forward_derivative, cur, k, dt) - f(cur));
                sum_3 += std::abs(maclaurin_series(f, back_derivative   , cur, k, dt) - f(cur));
            }
            ++count;
        }
        x.push_back(k);
        y_1.push_back(sum_1 / count);
        if (k != kmax) { 
            y_2.push_back(sum_2 / count);
            y_3.push_back(sum_3 / count);
        }
    }

    plot(x, y_1)->line_width(1).color("red");
    grid(true);
    show();
    plot(x, y_2)->line_width(1).color("blue");
    show();
    plot(x, y_3)->line_width(1).color("green");
    show();

    return 0;
}