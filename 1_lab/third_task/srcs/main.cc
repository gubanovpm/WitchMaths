#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <matplot/matplot.h>

using namespace boost::numeric::ublas;

#define EPSILON 1e-3
#define EPS 0

int main() {
    //========================================================================================
    // Табличные данные
    //========================================================================================
    matrix<double> z {3, 1}; 
    z(0, 0) = 0.0189, z(1, 0) = 0.1567, z(2, 0) = 0.8244;
    matrix<double> M  = {3, 1}; // кг/моль
    M(0, 0) = 0.028013, M(1, 0) = 0.04401, M(2, 0) = 0.03007;
    matrix<double> Pc = {3, 1}; // Па
    Pc(0, 0) = 3394400, Pc(1, 0) = 7386600, Pc(2, 0) = 4883900;
    matrix<double> Tc = {3, 1}; // К
    Tc(0, 0) = 126.2, Tc(1, 0) = 304.7, Tc(2, 0) = 305.43;

    matrix<double> OmegaA0 = {3, 1};
    OmegaA0(0, 0) = 0.45724, OmegaA0(1, 0) = 0.45724, OmegaA0(2, 0) = 0.45724; 
    matrix<double> OmegaB  = {3, 1};
    OmegaB(0, 0) = 0.077796, OmegaB(1, 0) = 0.077796, OmegaB(2, 0) = 0.077796;
    matrix<double> omega   = {3, 1};
    omega(0, 0) = 0.04, omega(1, 0) = 0.225, omega(2, 0) = 0.0986;

    matrix<double> k {3, 3};
    k(0, 0) = 0, k(0, 1) = -0.012, k(0, 2) = 0.1;
    k(1, 0) = -0.012, k(1, 1) = 0, k(1, 2) = 0.1;
    k(2, 0) = 0.1, k(2, 1) = 0.1, k(2, 2) = 0;

    matrix<double> lambdas = {3, 1};
    lambdas(0, 0) = 1, lambdas(1, 0) = 1, lambdas(1, 0) = 1;
    //========================================================================================
    // Пользовательские данные (вторая размерность =2 так как по заданию для двух пар (T, P))
    //========================================================================================
    matrix<double> T = {1, 1}; // К
    T(0, 0) = 250/*, T(1, 0) = 200, T(2, 0) = 200*/;
    

    matrix<double> P = {3, 1}; // Па
    P(0, 0) = 2000000, P(1, 0) = 2000000, P(2, 0) = 2000000;
    //==============================================================================================================================================
    // Вспомогательные лямбда-функции
    //==============================================================================================================================================
    auto get_K_i  = [](double omega_i, double T_ri, double P_ri) { return std::exp(5.372697*(3.6 + omega_i) * (1 - 1./T_ri)) / P_ri; };
    auto get_Tr_i = [](double T_i, double T_ci) { return T_i / T_ci; };
    auto get_Pr_i = [](double P_i, double P_ci) { return P_i / P_ci; };
    auto get_m_i  = [](double omega_i) { return 0.37964 + 1.408503 * omega_i - 0.16442*std::pow(omega_i, 2) + 0.016666 * std::pow(omega_i, 3); };
    auto get_OmegaA_i = [](double OmegaA0, double m_i, double T_rj) { return OmegaA0*std::pow(1+m_i*(1-std::sqrt(T_rj)), 2) ; };
    auto get_a_i  = [](double OmegaAi, double p_rj, double T_rj) { return OmegaAi * p_rj / std::pow(T_rj, 2) ; };
    auto get_b_i  = [](double OmegaBi, double p_rj, double T_rj) { return OmegaBi * p_rj / T_rj ; };
    auto get_a_ij = [](double k_ij, double a_i, double a_j) { return (1-k_ij) * std::sqrt(a_i * a_j); };
    auto get_b    = [](const matrix<double> &b, const matrix<double> &c) { 
        double sum = 0; 
        for (size_t i = 0; i < b.size1(); ++i)
            sum += b(i, 0) * c(i, 0);
        return sum; 
    };
    auto get_a   = [](const matrix<double> &a, const matrix<double> &c) {
        double sum = 0; 
        for (size_t i = 0; i < a.size1(); ++i) {
            for (size_t j = 0; j < a.size2(); ++j) {
                sum += a(i, j) * c(i, 0) * c(j, 0);
            }
        }
        return sum;
    };
    auto get_s_i = [](const matrix<double> &a, const matrix<double> &c, const size_t i) {
        double sum = 0;
        for (size_t k = 0; k < c.size1(); ++k) {
            sum += a(i, k) * c(k, 0);
        }
        return sum;
    };
    auto get_f_i = [](const double _p, const double _T, const double _c_i, const double _Z, const double _a, const double _b, const double _s_i, const double _a_i, const double _b_i) {
        return std::log(_p * _c_i) - std::log(_Z - _b) + _a / (2 * std::sqrt(2) * _b) * (2*_s_i/_a - _b_i/_b) * std::log((_Z - (std::sqrt(2) - 1) * _b) / (_Z + (std::sqrt(2) + 1) * _b)) + _b_i/_b * (_Z - 1);
    };

    auto function_alpha = [] (const matrix<double> &K, const matrix<double> &z, const double alpha) { 
        double sum = 0;
        for (size_t i = 0 ; i < K.size1(); ++i) {
            sum += ((z(i, 0) * (K(i, 0) - 1))/(alpha * (K(i, 0) - 1) + 1));
        }
        return sum;
    };

    //==============================================================================================================================================
    // Основная функция
    //==============================================================================================================================================
    std::vector<double> phi = {1.};
    // size_t count = 100;
    // for (size_t i = 0 ; i < count; ++i) phi.push_back(double(i)/count);

    std::vector<double> x = {};
    std::vector<size_t> y = {};
    
    for (size_t i = 0; i < phi.size(); ++i) {
        double PHI = phi[i];
        for (size_t pairs_count = 0 ; pairs_count < std::min(T.size1(), P.size1()); ++pairs_count) {
            // Сбросим значения для первой итерации
            z(0, 0) = 0.0189, z(1, 0) = 0.1567, z(2, 0) = 0.8244;
            size_t iteration = 0;
            matrix<double> c_g = {3, 1},  c_l = {3, 1};
            matrix<double> c_g_star = {3, 1}, c_l_star = {3, 1};
            matrix<double> K  = {3, 1};
            matrix<double> Tr = {3, 1};
            matrix<double> Pr = {3, 1};
            matrix<double> m  = {3, 1};
            matrix<double> OmegaA = {3, 1};
            matrix<double> a_i  = {3, 1};
            matrix<double> b_i  = {3, 1};
            
            matrix<double> a_ij = {3, 3};
            matrix<double> s_l_i = {3, 1}, s_g_i = {3, 1};
            
            // Вызов лямбда функций для заполнения соответствующих векторов
            for (size_t i = 0; i < 3; ++i) {
                Tr(i, 0) = get_Tr_i(T(pairs_count, 0), Tc(i, 0));
                Pr(i, 0) = get_Pr_i(P(pairs_count, 0), Pc(i, 0));
                m(i, 0) = get_m_i(omega(i, 0));
                OmegaA(i, 0) = get_OmegaA_i(OmegaA0(i, 0), m(i, 0), Tr(i, 0));
                K(i, 0)   = get_K_i(omega(i, 0), Tr(i, 0), Pr(i, 0));
                b_i(i, 0) = get_b_i(OmegaB(i, 0) , Pr(i, 0), Tr(i, 0));
                a_i(i, 0) = get_a_i(OmegaA0(i, 0), Pr(i, 0), Tr(i, 0));
            }
            // std::cout << K << std::endl;

            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    a_ij(i, j) = get_a_ij(k(i, j), a_i(i, 0), a_i(j, 0));
                }
            }

            while (true) {

                // рассчитаем коэффициент квадратоного уравнения (в данном случае это будут решения квадратного уравнения)
                double a_koef = (z(0, 0) + z(1, 0) + z(2, 0)) * (K(0, 0)-1) * (K(1, 0)-1) * (K(2, 0)-1) ;
                double b_koef = z(0, 0) * (K(0, 0) - 1) * (K(1, 0) + K(2, 0) - 2) +
                                z(1, 0) * (K(1, 0) - 1) * (K(0, 0) + K(2, 0) - 2) +
                                z(2, 0) * (K(2, 0) - 1) * (K(0, 0) + K(1, 0) - 2);
                double c_koef = z(0, 0) * (K(0, 0) - 1) + z(1, 0) * (K(1, 0) - 1) + z(2, 0) * (K(2, 0) - 1);
                double D = std::pow(b_koef, 2) - 4 * a_koef * c_koef;

                double xx1 = (-b_koef - std::sqrt(D)) / (2 * a_koef), xx2 = (-b_koef + std::sqrt(D)) / (2 * a_koef);

                // std::cout << "Корни по альфа: " <<  (-b_koef - std::sqrt(D)) / (2 * a_koef) << " " << (-b_koef + std::sqrt(D)) / (2 * a_koef) << std::endl;

                // Дихотомия
                double l_bound = 0, r_bound = 1, value , stop = 0;
                double l_sign = function_alpha(K, z, l_bound);
                double r_sign = function_alpha(K, z, r_bound);
                while (std::abs(value = function_alpha(K, z, (l_bound + r_bound) / 2)) > EPSILON && stop != 1000) {
                    if (value * l_sign > 0) l_bound = (l_bound + r_bound) / 2;
                    if (value * r_sign > 0) r_bound = (l_bound + r_bound) / 2;
                    ++stop;

                    // std::cout << (l_bound + r_bound) / 2 << std::endl;
                }

                std::cout << (l_bound + r_bound) / 2 << std::endl;
                double alpha = (l_bound + r_bound) / 2;
                //  double alpha = (std::abs(xx1) < 1 ? xx1 : xx2 );
                std::cout << "alpha = " << alpha << std::endl;
                // here's strange partition of code
                if (iteration == 0)
                for (size_t i = 0; i < 3; ++i) {
                    c_g(i, 0) = z(i, 0) * K(i, 0) / (alpha * (K(i,0) - 1) + 1);
                    c_l(i, 0) = z(i, 0) / (alpha * (K(i,0) - 1) + 1);
                }

                if (iteration != 0)
                for (size_t i = 0; i < 3; ++i) {
                    c_g_star(i, 0) = z(i, 0) * K(i, 0) / ( alpha * (K(i, 0) - 1) + 1);
                    c_g(i , 0) = PHI * c_g_star(i, 0) + (1 - PHI) * c_g(i , 0);

                    c_l_star(i, 0) = z(i, 0) / ( alpha * (K(i, 0) - 1) + 1);
                    c_l(i , 0) = PHI * c_l_star(i, 0) + (1 - PHI) * c_l(i, 0);
                }

                double a_l = get_a(a_ij, c_l), a_g = get_a(a_ij, c_g);
                double b_l = get_b(b_i , c_l), b_g = get_b(b_i , c_g);

                for (size_t i = 0; i < 3; ++i) {
                    s_l_i(i, 0) = get_s_i(a_ij, c_l, i);
                    s_g_i(i, 0) = get_s_i(a_ij, c_g, i);
                }

                // Решим кубическое уравнение тригонометрическим способом
                double a_c_l = b_l - 1, b_c_l = a_l - 2 * b_l - 3 * std::pow(b_l, 2), c_c_l = std::pow(b_l, 3) + std::pow(b_l, 2) - a_l * b_l;
                double Q_l = (3*b_c_l - std::pow(a_c_l, 2)) / 9;
                double R_l = (9*a_c_l*b_c_l - 2 * std::pow(a_c_l, 3) - 27 * c_c_l) / 54;
                double S_l = std::pow(Q_l,  3) - std::pow(R_l, 2);

                double a_c_g = b_g - 1, b_c_g = a_g - 2 * b_g - 3 * std::pow(b_g, 2), c_c_g = std::pow(b_g, 3) + std::pow(b_g, 2) - a_g * b_g;
                double Q_g = -(3*b_c_g - std::pow(a_c_g, 2)) / 9;
                double R_g = -(9*a_c_g*b_c_g - 2 * std::pow(a_c_g, 3) - 27 * c_c_g) / 54;
                double S_g = std::pow(Q_g,  3) - std::pow(R_g, 2);

                double psi_l = 0, psi_g = 0;
                double t_1 = 0, t_2 = 0;
                double x_1, x_2, x_3;
                std::vector<double> x = {};

                // std::cout << "x^3 + " <<  a_c_l << " x^2 + " << b_c_l << " x + " << c_c_l << std::endl;
                // std::cout << "x^3 + " <<  a_c_g << " x^2 + " << b_c_g << " x + " << c_c_g << std::endl; 

                if (S_l > EPS) {
                    psi_l = 1./3 * std::acos(R_l / std::sqrt(std::pow(Q_l, 3)));
                    // std::cout << "psi L " <<-2. * std::sqrt(Q_l) * std::cos(psi_l) << std::endl;
                    x_1 = -2. * std::sqrt(Q_l) * std::cos(psi_l) - a_c_l / 3.;
                    x_2 = -2. * std::sqrt(Q_l) * std::cos(psi_l + M_PI * 2./3) - a_c_l / 3.;
                    x_3 = -2. * std::sqrt(Q_l) * std::cos(psi_l - M_PI * 2./3) - a_c_l / 3.;
                    if (x_1 > 0) x.push_back(x_1);
                    if (x_2 > 0) x.push_back(x_2);
                    if (x_3 > 0) x.push_back(x_3);

                    // std::cout << x_1 << " ; " << x_2 << " ; " << x_3 << std::endl;
                } else if (S_l < -EPS) {
                    if (Q_l > EPS) {
                        psi_l = 1./3 * std::acosh(std::abs(R_l) / std::sqrt(std::abs(std::pow(Q_l, 3))));
                        double siln_R_l = (R_l == 0 ? 0. : R_l / std::abs(R_l));
                        x_2 = -2. * siln_R_l * std::sqrt(std::abs(Q_l)) * std::cosh(psi_l) - a_c_l / 3;

                        if (x_2 > 0) x.push_back(x_2);
                    }
                    else {
                        psi_l = 1./3 * std::asinh(std::abs(R_l) / std::sqrt(std::abs(std::pow(Q_l, 3))));
                        double siln_R_l = (R_l == 0 ? 0. : R_l / std::abs(R_l));
                        x_2 = -2. * siln_R_l * std::sqrt(std::abs(Q_l)) * std::sinh(psi_l) - a_c_l / 3;

                        if (x_2 > 0) x.push_back(x_2);
                    }                  
                }

                if (x.size() > 0)
                    t_1 = *std::min_element(x.begin(), x.end());
                x = {};
                if (S_g > EPS) {
                    psi_g = 1./3 * std::acos(R_g / std::sqrt(std::pow(Q_g, 3)));
                    x_1 = -2 * std::sqrt(Q_g) * std::cos(psi_g) - a_c_g/3;
                    x_2 = -2 * std::sqrt(Q_g) * std::cos(psi_g + M_PI * 2/3) - a_c_g / 3;
                    x_3 = -2 * std::sqrt(Q_g) * std::cos(psi_g - M_PI * 2/3) - a_c_g / 3;
                    if (x_1 > 0) x.push_back(x_1);
                    if (x_2 > 0) x.push_back(x_2);
                    if (x_3 > 0) x.push_back(x_3);

                } else if (S_g < -EPS) {
                    if (Q_g > EPS) {
                        psi_g = 1./3 * std::acosh(std::abs(R_g) / std::sqrt(std::abs(std::pow(Q_g, 3))));
                        double sign_R_g = (R_g == 0 ? 0. : R_g / std::abs(R_g));
                        x_2 = -2. * sign_R_g * std::sqrt(std::abs(Q_g)) * std::cosh(psi_g) - a_c_g / 3;

                        if (x_2 > 0) x.push_back(x_2);
                    }
                    else {
                        psi_g = 1./3 * std::asinh(std::abs(R_g) / std::sqrt(std::abs(std::pow(Q_g, 3))));
                        double sign_R_g = (R_g == 0 ? 0. : R_g / std::abs(R_g));
                        x_2 = -2. * sign_R_g * std::sqrt(std::abs(Q_g)) * std::sinh(psi_g) - a_c_g / 3;

                        if (x_2 > 0) x.push_back(x_2);
                    }
                }

                if (x.size() > 0)
                    t_2 = *std::max_element(x.begin(), x.end());


                std::cout << t_1 << " ; " << t_2 << std::endl;
                // for (size_t i = 0; i < 3 ; ++i)
                //     std::cout << c_g(i, 0) << std::endl;
                // for (size_t i = 0; i < 3 ; ++i)
                //     std::cout << c_l(i, 0) << std::endl;

                // std::cout << Q_l << " ; " << R_l << " ; "  << S_l << " >> " << psi_l << " " << t_1 << std::endl;
                // std::cout << Q_g << " ; " << R_g << " ; "  << S_g << " >> " << psi_g << " " << t_2 << std::endl;

                matrix<double> f_l = {3, 1}, f_g = {3, 1};
                for (size_t i = 0; i < 3; ++i) {
                    f_l(i, 0) = get_f_i(P(pairs_count, 0), T(pairs_count, 0), c_l(i, 0), t_1, a_l, b_l, s_l_i(i, 0), a_i(i, 0), b_i(i, 0));
                    f_g(i, 0) = get_f_i(P(pairs_count, 0), T(pairs_count, 0), c_g(i, 0), t_2, a_g, b_g, s_g_i(i, 0), a_i(i, 0), b_i(i, 0));
                    K(i, 0) *= (f_l(i, 0) / f_g(i, 0));
                }

                std::cout << c_l << std::endl << c_g << std::endl << z << std::endl << std::endl;;
                bool flag = true;
                for (size_t i = 0; i < 3; ++i) {
                    flag &= (std::abs(f_l(i, 0) / f_g(i, 0) - 1) < EPSILON);
                    if (!flag) break;
                }

                if (flag || (iteration == 20)) break;
                ++iteration;

                // if (iteration == 3) break;
            }
            std::cout << "liq : " << c_l << std::endl;
            std::cout << "gaze : " << c_g << std::endl << std::endl;
            x.push_back(PHI);
            y.push_back(iteration);
        }
    }

    // matplot::plot(x, y)->line_width(1).color("red");
    // matplot::xlabel("Параметр Фи");
    // matplot::ylabel("Количество итераций");
    // matplot::show();
    
    return 0;
}