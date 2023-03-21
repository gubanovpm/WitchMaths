#include "../libs/runge_kutta.hh"

/** Функция правой части дифференциального уравнения y' = f(t, y)
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */ 
Vec rightPart(const Time& time, const Vec& state) noexcept {
    Vec result(2);
    // std::cout << "I am in function" << std::endl;
    result(0) = state(1);
    result(1) = std::cos(3*time) - 4 * state(0);
    // std::cout << "Getting outfrom function" << std::endl;
    return result;
}

int main() {
    ButcherTable<4> table;
    table.column = {   0, 1./2, 1./2,    1};
    table.string = {1./6, 1./3, 1./3, 1./6};
    table.matrix = std::array<std::array<double, 4>, 4>  {
        std::array<double, 4>{  0,   0, 0, 0}, 
        std::array<double, 4>{0.5,   0, 0, 0}, 
        std::array<double, 4>{  0, 0.5, 0, 0}, 
        std::array<double, 4>{  0,   0, 1, 0}
    };
    double step = 0.01;
    State state(Vec(2), 0); state.state << 2, -2.2;
    unsigned iterations = 100;

    std::vector<Vec> runge_kutta = explicitRK(state, step, iterations, rightPart, table);
    for (unsigned i = 0; i < runge_kutta.size(); ++i) {
        std::cout << runge_kutta[i] << std::endl;
    }

    return 0;
}