#include "seidel.hh"

int main () {
    size_t n_upper      = 10 ;
    size_t seidel_upper = 1000;

    double dt = 0.001;
    double dx = 0.001;

    std::vector<size_t> n_vector = {};
    std::vector<double> a_vector = {};
    std::vector<size_t> k_vector = {};

    for (size_t n = 1; n <= n_upper; ++n) {

        for (double a = 0; a <= 1; a += dt) {
            matrix<double> A   = get_right_matrix(n, a);
            matrix<double> f   = get_right_vector(n, a);
            matrix<double> ans = get_right_answer(n, a);

            seidel_method seidel(A, f);

            size_t itter_count = 0;
            double delta = dx * 2;
            do {
                if (itter_count == seidel_upper) break;
                delta = std::abs(norm().euclid(ans) - norm().euclid(seidel.get_k_itter(itter_count)));
                ++itter_count; 
            } while (delta > dx); --itter_count;

            if (itter_count < 1000) {
                n_vector.push_back(n);
                a_vector.push_back(a);
                k_vector.push_back(itter_count);
            }
        }
    }

    auto l = matplot::plot3(n_vector, a_vector, k_vector)->line_width(1).color("red");
    matplot::xlabel("Размер системы уравнений");
    matplot::ylabel("Значение параметра a");
    matplot::zlabel("Количество необходимых иттераций метода Зейделя для погрешности = 0.001");
    matplot::show();

    return 0;
}