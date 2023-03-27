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
                unsigned count_of_it = 15;
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
        if (isImplicit(table)) return implicitRK(initial, step, iterations, rightPart, table);
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
 * @param a  Левая граница отрезка
 * @param b  Правая граница отрезка
 * @param y0 y(a)
 * @param y1 y(b)
 * @param h Шаг интегрирования
 * @param eps Точность
 * @param table Таблица Бутчера для вызова метода Рунге-Кутты
 * 
 * @return массив решений
 */
// template <unsigned s>
// std::vector<Vec> shooting_method(
//     const double a, 
//     const double b, 
//     const double y0, 
//     const double y1, 
//     const double h, 
//     const double eps,
//     const ButcherTable<s>& table
//     ) noexcept {
//         double nu1 = 1.0, nu2 = 0.8, nu3; // некоторое значение тангенса угла наклона касательной к решению в точке x = a
//         int N = (int)((b - a) / h + 1); // количество элементов
//         double[,] Y1 = new double[2, N], Y2 = new double[2, N]; // массивы для хранения решений
//         Y1 = rk4(a, b, h, y0, nu1, Y1); // решаем систему для nu1
//         Y2 = rk4(a, b, h, y0, nu2, Y2); // решаем систему для nu2
//         while (Math.Abs(Y2[0, N - 1] - y1) > eps) //
//         {
//             nu3 = nu2 - ((nu2 - nu1) / (Y2[0, N - 1] - Y1[0, N - 1])) * (Y2[0, N - 1] - y1); // вычисляем новое nu
//             nu1 = nu2;
//             nu2 = nu3;
//             Y1 = Y2;
//             Y2 = new double[2, N]; // без этого не работает ?!? о_О
//             Y2 = rk4(a, b, h, y0, nu2, Y2); //вычисляем решение для нового nu
//         }
//         return Y2;
//     }
 

#endif