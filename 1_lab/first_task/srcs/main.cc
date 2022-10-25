#include "maclaurin.hh"

using namespace matplot;

#define EXP

int main() {
    #ifdef SIN
        auto f = [](double x){ return std::sin(x); };
    #else 
        auto f = [](double x){ return std::exp(x); };
    #endif
    
    maclaurin_series msc(f, maclaurin_series::CENTRAL);
    maclaurin_series msb(f, maclaurin_series::BACK   );
    maclaurin_series msf(f, maclaurin_series::FORWARD);

    double xbegin = 0, xend = 1, dt = 0.001;
    size_t kbegin = 1, dk = 1;

    std::vector<maclaurin_series> ms = {};
    ms.push_back(msc);
    ms.push_back(msb);
    ms.push_back(msf);
    std::vector<std::string> mttl = {"Central", "Back", "Forward"};
    std::vector<size_t> akmax = {7, 6, 6};


    for (size_t i = 0; i < ms.size(); ++i) {
        std::vector<size_t> x = {};
        std::vector<double> y = {};

        
        for (size_t k = kbegin; k <= akmax[i]; k += dk) {
            double sum = 0; size_t count = 0;

            for (double t = xbegin; t <= xend; t += dt) {
            // std::cout << "I AM GAY" << std::endl;
                sum += std::abs(ms[i].value(k, t, dt) - f(t));
                ++count;
            }
            x.push_back(k);
            y.push_back(sum / count);
        }

        title(mttl[i]);
        grid(true);
        xlabel("Количество слагаемых в ряду Маклорена n, шт");
        ylabel("Среднее тклонение по норме abs");
        plot(x, y)->line_width(1).color("red");
        show();
    }


    return 0;
}