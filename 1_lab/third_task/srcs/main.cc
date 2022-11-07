#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

#define EPSILON 1e-6

int main() {
    //========================================================================================
    // Табличные данные
    //========================================================================================
    matrix<double> z {3, 1}; 
    z(0, 0) = 0.0189, z(1, 0) = 0.1567, z(2, 0) = 0.8244;
    matrix<double> M  = {3, 1}; // кг/моль
    M(0, 0) = 0.028013, M(1, 0) = 0.04401, M(2, 0) = 0.03007;
    matrix<double> Pc = {3, 1}; // Па
    Pc(0, 0) = 0.00033944, Pc(1, 0) = 0.00073866, Pc(2, 0) = 0.00048839;
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
    matrix<double> T = {3, 1}; // К
    T(0, 0) = 73, T(1, 0) = 73, T(2, 0) = 73;

    matrix<double> P = {3, 1}; // Па
    P(0, 0) = 0.0001, P(1, 0) = 0.0001, P(2, 0) = 0.0001;
    //==============================================================================================================================================
    // Вспомогательные лямбда-функции
    //==============================================================================================================================================
    auto get_K_i  = [](double omega_i, double T_ri, double P_ci) { return std::exp(5.373*(1-omega_i)) * (1 - 1/T_ri) * P_ci; };
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
    auto get_s_i = [](const matrix<double> &a, const matrix<double> &c, const int i) {
        double sum = 0;
        for (size_t k = 0; k < c.size1(); ++k) {
            sum += a(i, k) * c(k, 0);
        }
        return sum;
    };
    auto get_f_i = [](const double _p, const double _T, const double _c_i, const double _Z, const double _a, const double _b, const double _s_i, const double _a_i, const double _b_i) {
        return std::log(_p * _c_i) + std::log(_Z - _b) + _a / (2 * std::sqrt(2) * _b) * (2*_s_i/_a - _b_i/_b) * std::log((_Z - (std::sqrt(2) - 1) * _b) / (_Z + (std::sqrt(2) + 1) * _b)) + _b_i/_b * (_Z - 1);
    };

    //==============================================================================================================================================
    // Основная функция
    //==============================================================================================================================================
    double PHI = 0;
    for (size_t pairs_count = 0; pairs_count < std::min(T.size1(), P.size1()); ++pairs_count) {
        // Сбросим значения для первой итерации
        z(0, 0) = 0.0189, z(1, 0) = 0.1567, z(2, 0) = 0.8244;
        size_t iteration = 0;
        matrix<double> c_g = {3, 1}, c_l = {3, 1};
        matrix<double> c_g_star = {3, 1}, c_l_star = {3, 1};
        matrix<double> K  = {3, 1};
        matrix<double> Tr = {3, 1};
        matrix<double> Pr = {3, 1};
        matrix<double> m  = {3, 1};
        matrix<double> OmegaA = {3, 1};
        matrix<double> a_i  = {3, 1};
        matrix<double> b_i  = {3, 1};
        
        matrix<double> a_ij = {3, 3};
        matrix<double> s_i = {3, 1};
        
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
            // std::cout << "Корни по альфа: " <<  (-b_koef - std::sqrt(D)) / (2 * a_koef) << " " << (-b_koef + std::sqrt(D)) / (2 * a_koef) << std::endl;

            // I don't know how we are must calc two roots of equation
            double Z_l = 1, Z_g = 1;

            double alpha = (-b_koef + std::sqrt(D)) / (2 * a_koef);
            // here's strange partition of code
            c_g(0, 0) = z(0, 0) * K(0, 0) / (alpha * (K(0,0) - 1) + 1);
            c_l(0, 0) = z(0, 0) / (alpha * (K(0,0) - 1) + 1);
            
            for (size_t i = 1; i < 3; ++i) {
                c_g_star(i, 0) = z(i, 0) * K(i, 0) / ( lambdas(i, 0) * (K(i, 0) - 1) + 1);
                c_g(i , 0) = PHI * c_g_star(i, 0) + (1 - PHI) * c_g(i - 1, 0);

                c_l_star(i, 0) = z(i, 0) / ( lambdas(i, 0) * (K(i, 0) - 1) + 1);
                c_l(i , 0) = PHI * c_l_star(i, 0) + (1 - PHI) * c_l(i - 1, 0);
            }

            double a_l = get_a(a_ij, c_l), a_g = get_a(a_ij, c_g);
            double b_l = get_b(b_i , c_l), b_g = get_b(b_i , c_g);

            matrix<double> f_l = {3, 1}, f_g = {3, 1};
            for (size_t i = 0; i < 3; ++i) {
                f_l(i, 0) = get_f_i(P(pairs_count, 0), T(pairs_count, 0), c_l(i, 0), Z_l, a_l, b_l, s_i(i, 0), a_i(i, 0), b_i(i, 0));
                f_g(i, 0) = get_f_i(P(pairs_count, 0), T(pairs_count, 0), c_g(i, 0), Z_g, a_g, b_g, s_i(i, 0), a_i(i, 0), b_i(i, 0));
                K(i, 0) *= (f_l(i, 0) / f_g(i, 0));
            }

            bool flag = true;
            for (size_t i = 0; i < 3; ++i) {
                flag = (std::abs(f_l(i, 0) / f_g(i, 0) - 1) < EPSILON);
                if (!flag) break;
            }

            if (flag || (iteration == 10000)) break;
            ++iteration;
        }
        std::cout << iteration << std::endl;
    }
    
    return 0;
}