#ifndef __diff_eq_hh__
#define __diff_eq_hh__

#include <array>
#include <cmath>
#include <vector>
#include <functional>
#include <eigen3/Eigen/Core>


#include <sciplot/sciplot.hpp>
#include <iostream>

using Vec  = Eigen::VectorXd;
using Mat  = Eigen::MatrixXd;
using Time = double;

struct State {
    Vec state;
    Time t;
};

/** Таблица Бутчера **/
template <unsigned s>
struct ButcherTable {
    std::array<double, s> column;
    std::array<double, s> string;
    std::array<std::array<double, s>, s> matrix;
};

/** Функция явного метода Рунге-Кутты
 * @tparam s  стадийность метода
 * @param initial начальное условие
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param table таблица Бутчера метода
 * 
 * @return массив решений
 */
template <unsigned s>
std::vector<Vec> explicitRK(
    const State &initial,
    double step,
    unsigned iterations,
    const std::function<Vec(Time, Vec)> &rightPart,
    const ButcherTable<s>& table
    ) noexcept {
        std::vector<Vec> result;
        result.push_back(initial.state);

        for (unsigned i = 1; i <= iterations; ++i) {
            std::vector<Vec> k = {};
            for (unsigned s_i = 0; s_i < s; ++s_i) {
                double x_n = initial.t + step * table.column[s_i] + step * (i - 1);
                Vec y_n = result[i-1];
                for (unsigned s_j = 0; s_j < s_i; ++s_j) {
                    y_n += (step * table.matrix[s_i][s_j] * k[s_j]); 
                }
                k.push_back(rightPart(x_n, y_n));
            }
            
            Vec temp = result[i-1];
            for (unsigned s_i = 0; s_i < s; ++s_i) {
                temp += step * table.string[s_i] * k[s_i];
            }
            result.push_back(temp);
        }

        return result;
    }

/** Функция неявного метода Рунге-Кутты
 * @tparam s  стадийность метода
 * @param initial начальное условие
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param table таблица Бутчера метода
 *
 * @return массив решений
 */
template <unsigned s>
std::vector<Vec> implicitRK(
    const State& initial,
    double step,
    unsigned iterations,
    const std::function<Vec(Time, Vec)>& rightPart,
    const ButcherTable<s>& table
    ) noexcept {
        std::vector<Vec> result;
        result.push_back(initial.state);
        for (unsigned i = 1; i <= iterations; ++i) {
            std::vector<Vec> k = {};
            for (unsigned j = 0; j < s; ++j) k.push_back(Vec::Zero(initial.state.size()));

            for (unsigned s_i = 0; s_i < s; ++s_i) {
                double x_n = initial.t + step * table.column[s_i] + step * (i - 1);

                // Метод простой иттерации для коэффициентов Рунге-Кутты
                // Количество иттераций для простого метода иттераций
                unsigned count_of_it = 20;
                for (unsigned i_it = 1; i_it <= count_of_it; ++i_it) {
                    for (unsigned s_i_it = 0; s_i_it < s; ++s_i_it) {
                        double t_n = initial.t + step * table.column[s_i_it] + step * (i_it - 1);
                        Vec tmp = result[i-1];
                        for (unsigned s_j_it = 0; s_j_it < s; ++s_j_it) {
                            tmp += (step * table.matrix[s_i_it][s_j_it] * k[s_j_it]);
                        }
                        k[s_i_it] = rightPart(t_n, tmp);
                    }
                }
                // Конец метода простых иттераций

                Vec y_n = result[i-1];
                for (unsigned s_j = 0; s_j < s; ++s_j) {
                    y_n += (step * table.matrix[s_i][s_j] * k[s_j]); 
                }
                k[s_i] = rightPart(x_n, y_n);
            }
            
            Vec temp = result[i-1];
            for (unsigned s_i = 0; s_i < s; ++s_i) {
                temp += step * table.string[s_i] * k[s_i];
            }
            result.push_back(temp);
        }

        return result;
    }
/** Функция неявного метода Адамса
 * @tparam s порядок метода
 * @param table Таблица Бутчера для вызова метода Рунге-Кутты
 * 
 * @return массив решений
 */
template <unsigned s>
bool isImplicit(const ButcherTable<s>& table) {
    for (unsigned i = 0; i < s; ++i) {
        for (unsigned j = i; j < s; ++j) {
            if (table.matrix[i][j]) {
                return true;
            }
        }
    }
    return false;
}

/** Функция явного метода Рунге-Кутты
 * @tparam s  стадийность метода
 * @param initial начальное условие
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param table таблица Бутчера метода
 * 
 * @return массив решений
 */
template <unsigned s>
std::vector<Vec> RungeKutta(
    const State& initial,
    double step,
    unsigned iterations,
    const std::function<Vec(Time, Vec)>& rightPart,
    const ButcherTable<s>& table
    ) noexcept {
        if (isImplicit(table)) { 
            return implicitRK(initial, step, iterations, rightPart, table);
        }
        return explicitRK(initial, step, iterations, rightPart, table);
    }

/** Функция явного метода Адамса
 * @tparam s порядок метода
 * @param initial начальное условие (после разгона)
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param coefs коэффициенты при правых частях
 * @param previousRightParts - правые части дифференциального уравнения до начального условия (разгон метода)
 *
 * @return массив решений
 */
template <unsigned s>
std::vector<Vec> explicitAdams(
    const State& initial,
    double step,
    unsigned iterations,
    const std::function<Vec(Time, Vec)>& rightPart,
    const std::array<double, s>& coefs,
    const std::array<Vec, s>& previousRightParts
    ) noexcept {
        std::vector<Vec> result;
        for (unsigned i = 0; i < s; ++i)
            result.push_back(coefs[i]);
        for (unsigned i = s; i <= iterations; ++i) {
            Vec temp = result[i - 1];
            for (unsigned j = 1; j <= s; ++j) {
                double x_n = initial.t + step * table.column[s_i] + step * (i - j - 1);
                temp += step * coefs[j] * rightPart(x_n, result[i - j - 1]);
            }
            result.push_back(temp);
        }
        return result;
    }

/** Функция неявного метода Адамса
 * @tparam s порядок метода
 * @param initial начальное условие (после разгона)
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param coefs коэффициенты при правых частях
 * @param previousRightParts - правые части дифференциального уравнения до начального условия (разгон метода)
 *
 * @return массив решений
 */
// template <unsigned s>
// std::vector<Vec> implicitAdams(
//     const State& initial,
//     double step,
//     unsigned iterations,
//     const std::function<Vec(Time, Vec)>& rightPart,
//     const std::array<double, s + 1>& coefs,
//     const std::array<Vec, s>& rightParts
//     ) noexcept {

//     }

/** Функция неявного метода Адамса
 * @tparam s порядок метода
 * @param segment границы отрезка
 * @param state значения функции на концах отрезка
 * @param tangens начальные значения угла наклона
 * @param rightPart функция правой части
 * @param table таблица Бутчера для вызова метода Рунге-Кутты
 * @param h шаг интегрирования
 * @param eps точность
 * 
 * @return массив решений
 */
template <unsigned s>
std::vector<Vec> shooting_method(
    const Vec &segment,
    const Vec &state,
    const Vec &tangens,
    const std::function<Vec(Time, Vec)>& rightPart,
    const ButcherTable<s>& table,
    const double h,
    const double eps
    ) noexcept {
        if (segment.size() != 2 || state.size() != 2 || tangens.size() != 2) {
            std::cout << "Wrong begin state size" << std::endl;
            return {};
        }

        int iterations = ((segment(1) - segment(0)) / h);
        double tg = tangens(1);

        State st_1; st_1.state = Vec(2); 
        st_1.state << state(0), tangens(0);
        st_1.t     = segment(0);
        
        State st_2; st_2.state = Vec(2); 
        st_2.state << state(0), tangens(1);
        st_2.t     = segment(0);
        
        State st_3; st_3.state = Vec(2); st_2.t = segment(0);

        std::vector<Vec> runge_kutta_1 = RungeKutta(st_1, h, iterations, rightPart, table);
        std::vector<Vec> runge_kutta_2 = RungeKutta(st_2, h, iterations, rightPart, table);
        
        unsigned it_count = 0;
        while (std::abs(runge_kutta_2[iterations](0) - state(1)) > eps) {
            tg = st_2.state(1) - (st_2.state(1) - st_1.state(1)) / (runge_kutta_2[iterations](0) - runge_kutta_1[iterations](0)) * (runge_kutta_2[iterations ](0) - state(1));
            st_3.state << state(0), tg ;
            st_1.state << state(0), st_2.state(1) ;
            st_2.state << state(0), tg ;
            runge_kutta_1 = runge_kutta_2;

            runge_kutta_2 = RungeKutta(st_2, h, iterations, rightPart, table);
            if (it_count == 1000) break;
            ++it_count;
        }
        std::cout << "Угловой коэффициент для начальной точки: " << tg << std::endl;
        std::cout << "Номер иттерации: " << it_count << std::endl;
        return runge_kutta_2;
    }
 

#endif