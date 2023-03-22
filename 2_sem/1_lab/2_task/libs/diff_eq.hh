#ifndef __diff_eq_hh__
#define __diff_eq_hh__

#include <array>
#include <vector>
#include <functional>
#include <eigen3/Eigen/Core>


#include <sciplot/sciplot.hpp>
#include <iostream>

using Vec  = Eigen::VectorXd;
using Mat  = Eigen::MatrixXd;
using Time = double;

using FunctionT = std::function<Vec(const Time &, const Vec &)>;

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

/** Функция явного метода рунге-Кутты
 * @tparam s  стадийность метода
 * @param initial начальное условие
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param table таблица Бутчера метода
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

/** Функция неявного метода рунге-Кутты
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
                        // for (unsigned j = 0; j < k.size(); ++j) {
                        //     std::cout << k[j].transpose() << std::endl;
                        // }
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

/** Функция явного метода адамса
 * @tparam s порядок метода
 * @param initial начальное условие (после разгона)
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param coefs коэффициенты при правых частях
 * @param previousRightParts - правые части дифференциального уравнения до начального условия (разгон метода)
 *
 * @return массив решений
 *
 * Методы Адамса требуют разгона. То есть необходимо несколько решений найти при помощи одношагового метода, а потом
 * только запустить многошаговый
 */
// template <unsigned s>
// std::vector<Vec> explicitAdams(
//     const State& initial,
//     double step,
//     unsigned iterations,
//     const std::function<Vec(Time, Vec)>& rightPart,
//     const std::array<double, s>& coefs,
//     const std::array<Vec, s>& previousRightParts
//     ) noexcept {

//     }

/** Функция неявного метода адамса
 * @tparam s порядок метода
 * @param initial начальное условие (после разгона)
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param coefs коэффициенты при правых частях
 * @param previousRightParts - правые части дифференциального уравнения до начального условия (разгон метода)
 *
 * @return массив решений
 *
 * Методы Адамса требуют разгона. То есть необходимо несколько решений найти при помощи одношагового метода, а потом
 * только запустить многошаговый
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

#endif