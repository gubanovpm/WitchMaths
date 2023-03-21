#ifndef __runge_kutta_hh__
#define __runge_kutta_hh__

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

    State(
        const Vec &s, 
        const Time t
    ):
        state(s), 
        t(t) {
        }
};

/** Таблица Бутчера **/
template <unsigned s>
struct ButcherTable {
    std::array<double, s> column;
    std::array<double, s> string;
    std::array<std::array<double, s>, s> matrix;
};

/** Функция явного метода рунге-Кутты
 * @tparam S  стадийность метода
 * @param initial начальное условие
 * @param step  шаг интегрирования
 * @param iterations количество шагов, которые необходимо сделать
 * @param rightPart функция правой части дифференциального уравнения
 * @param table таблица бутчера метода
 * @return массив решений
 */
template <unsigned S>
std::vector<Vec> explicitRK(
    const State &initial,
    double step,
    unsigned iterations,
    const std::function<Vec(Time, Vec)> &rightPart,
    const ButcherTable<S>& table
    ) noexcept {
        std::vector<Vec> result;
        result.push_back(initial.state);
        for (unsigned i = 1; i <= iterations; ++i) {
            std::vector<Vec> K = {};
            for (unsigned s_i = 0; s_i < S; ++s_i) {
                double t_n = initial.t + step * s_i;
                Vec y_n = Vec::Zero(initial.state.size());
                // std::cout << "Initial state size = " << initial.state.size() << std::endl;
                for (unsigned s_j = 0; s_j < s_i; ++s_j) {
                    Vec tmp = (step * table.matrix[s_i][s_j] * K[s_j]); 
                    // std::cout << "Current temp = " << tmp << std::endl;
                    // std::cout << table.matrix[s_i][s_j] << " " ;
                    for (unsigned j = 0; j < initial.state.size(); ++j) {
                        y_n += tmp;
                    }
                }
                // std::cout << "Current state y_n is : " << y_n.transpose() << std::endl;
                K.push_back(rightPart(t_n, y_n));
            }
            // for (unsigned j = 0; j < K.size(); ++j)
            //     std::cout << j << ":" << std::endl << K[i] << std::endl;
            Vec temp(initial.state.size());
            temp = result[i-1];
            for (unsigned s_i = 0; s_i < S; ++s_i) {
                temp += step * table.string[s_i] * K[s_i];
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
 * @param table таблица бутчера метода
 *
 * @return массив решений
 */
// template <unsigned s>
// std::vector<Vec> implicitRK(
//     const State& initial,
//     double step,
//     unsigned iterations,
//     const std::function<Vec(Time, Vec)>& rightPart,
//     const ButcherTable<s>& table
//     ) noexcept {
//         std::vector<Vec> result;
//         result.push_back(initial.state);
//         Vec K(S);

//         for (unsigned i = 1; i <= iterations; ++i) {
//             for (unsigned s_i = 0; s_i < S; ++s_i) {
//                 double t_n = initial.time + step * s_i;
//                 Vec y_n = result[s_i];
//                 for (unsigned s_j = 0; s_j < S; ++s_j) {
//                     y_n += h * table.matrix[s_i][s_j] * K(s_j);
//                 }
//                 K(s_i) = rightPart(t_n, y_n);
//             }
//             Vec temp(initial.state.size()) = result[i-1];
//             for (unsigned s_i = 0; s_i < S; ++s_i) {
//                 temp(s_i) += step * table.string[s_i] * K(s_i);
//             }
//             result.push_back(temp);
//         }

//         return result;
//     }

// /** Функция явного метода адамса
//  * @tparam s порядок метода
//  * @param initial начальное условие (после разгона)
//  * @param step  шаг интегрирования
//  * @param iterations количество шагов, которые необходимо сделать
//  * @param rightPart функция правой части дифференциального уравнения
//  * @param coefs коэффициенты при правых частях
//  * @param previousRightParts - правые части дифференциального уравнения до начального условия (разгон метода)
//  *
//  * @return массив решений
//  *
//  * Методы Адамса требуют разгона. То есть необходимо несколько решений найти при помощи одношагового метода, а потом
//  * только запустить многошаговый
//  */
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

// /** Функция неявного метода адамса
//  * @tparam s порядок метода
//  * @param initial начальное условие (после разгона)
//  * @param step  шаг интегрирования
//  * @param iterations количество шагов, которые необходимо сделать
//  * @param rightPart функция правой части дифференциального уравнения
//  * @param coefs коэффициенты при правых частях
//  * @param previousRightParts - правые части дифференциального уравнения до начального условия (разгон метода)
//  *
//  * @return массив решений
//  *
//  * Методы Адамса требуют разгона. То есть необходимо несколько решений найти при помощи одношагового метода, а потом
//  * только запустить многошаговый
//  */
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

// /** Функция неявного дифференцирования назад
//  * @tparam s порядок метода
//  * @param initial начальное условие (после разгона)
//  * @param step  шаг интегрирования
//  * @param iterations количество шагов, которые необходимо сделать
//  * @param rightPart функция правой части дифференциального уравнения
//  * @param coefs коэффициенты при решениях на текущем и предыдущих шагах
//  * @param previousRightParts - решения дифференциального уравнения до начального условия (разгон метода)
//  *
//  * @return массив решений
//  *
//  * Методы Адамса требуют разгона. То есть необходимо несколько решений найти при помощи одношагового метода, а потом
//  * только запустить многошаговый
//  */
// template <unsigned s>
// std::vector<Vec> implBackDif(
//     const State& initial,
//     double step,
//     unsigned iterations,
//     const std::function<Vec(Time, Vec)>& rightPart,
//     const std::array<double, s + 1>& coefs,
//     const std::array<Vec, s>& solutions
//     ) noexcept {

//     }


// /** Функция вложенных методов
//  * @tparam s стадийность
//  * @param initial начальное условие (после разгона)
//  * @param intialStep начальный шаг интегрирования
//  * @param maxTime максимальное время, до которого нужно интегрировать
//  * @param rightPart функция правой части дифференциального уравнения
//  * @param EmbeddedTable таблица Бутчера для вложенных методов
//  * @param errorFunction - функция оценки ошибки (оцениват ошибку по двум состояниям)
//  *
//  * @return массив решений
//  */
// template <unsigned s>
// std::vector<Vec> embeddedMethod(
//     const State& initial,
//     double intialStep,
//     Time maxTime,
//     const std::function<Vec(Time, Vec)>& rightPart,
//     const EmbeddedTable<s>& table,
//     double tolerance,
//     const std::function<double(Vec,Vec)>& errorFunction
//     ) noexcept {

//     }

#endif