#include "../libs/runge_kutta.hh"

/** Функция правой части дифференциального уравнения y' = f(t, y)
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */ 
Vec rightPart(const Time& time, const Vec& state) noexcept {
    Vec result(state.size());

    result(0) = state(1);
    result(1) = std::cos(3*time) - 4 * state(0);

    return result;
}

int main() {
    ButcherTable<4> table;
    table.column = {   0, 1./2, 1./2,    1};
    table.string = {1./6, 1./3, 1./3, 1./6};
    std::cout << std::endl;
    table.matrix = std::array<std::array<double, 4>, 4>  {
        std::array<double, 4> {    0,    0,    0,    0}, 
        std::array<double, 4> { 1./2,    0,    0,    0}, 
        std::array<double, 4> {    0, 1./2,    0,    0}, 
        std::array<double, 4> {    0,    0,    1,    0}
    };

    double step = 0.1;
    State state(Vec(2), 0); state.state << 0.8, 2;
    unsigned iterations = 80;

    std::vector<Vec> runge_kutta = explicitRK(state, step, iterations, rightPart, table);

    sciplot::Plot2D plot;
    plot.legend().atOutsideBottom().displayHorizontal().displayExpandWidthBy(2);

    sciplot::Vec x   = sciplot::linspace(0, step*iterations, iterations);
    sciplot::Vec y_1 = sciplot::linspace(0, 0, iterations);
    sciplot::Vec y_2 = sciplot::linspace(0, 0, iterations);
    for (unsigned i = 0; i < runge_kutta.size(); ++i) {
        y_1[i] = runge_kutta[i](0);
    }
    for (unsigned i = 0; i < runge_kutta.size(); ++i) {
        y_2[i] = runge_kutta[i](1);
    }
    plot.drawCurve(x, y_1).label("z(x)").lineWidth(2);
    plot.drawCurve(x, y_2).label("dz/dx(x)").lineWidth(2);
    plot.grid().lineWidth(2).show();
    sciplot::Figure figure = {{ plot }};
    sciplot::Canvas canvas = {{ figure }};
    canvas.size(720, 400);
    canvas.save("runge_kutta.png");

    return 0;
}