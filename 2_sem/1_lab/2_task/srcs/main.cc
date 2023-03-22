#include "../libs/diff_eq.hh"

double g(const double t) {
    return 9.81 + 0.01 * std::cos(2 * M_PI * t);
}

double L = 4;
double m = 2;
/** Функция правой части дифференциального уравнения y' = f(t, y)
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */ 
Vec rightPart(const Time& time, const Vec& state) noexcept {
    Vec result(state.size());
// Значение функций опредеяем здесь:
    // x y u v T
    // 0 1 2 3 4 
    result(0) = -state(1) / state(0) * state(3);
    result(1) = state(0) / L * std::sqrt(state(4) * L / m + state(1) * g(time));
    result(2) = -state(0) / m / L / state(4) ;
    result(3) =  -state(1) / m / L / state(4) - g(time);
    result(4) = -m*0.01*2*M_PI*std::sin(2*M_PI*time) * state(0) / L - m * g(time) * state(1) * state(3) / state(0) / L;

    return result;
}

int main() {
    // Таблица Бутчера для явного метода Рунге-Кутты 4 порядка
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

    // Таблица Бутчера для неявного метода Рунге-Кутты 3 порядка
    // ButcherTable<3> table;
    // table.column = {      0.32,  0.962963,  0.962963};
    // table.string = {  0.720046,  0.720046,  0.008391};
    // std::cout << std::endl;
    // table.matrix = std::array<std::array<double, 3>, 3>  {
    //     std::array<double, 3> {  0.333333,         0, -0.013333}, 
    //     std::array<double, 3> {  0.625153,  0.333333,  0.004477}, 
    //     std::array<double, 3> {  9.516331, -8.886702,  0.333333} 
    // };

    // Проинициализируем начальные значения
    State state;
    state.state = Vec(5); state.state << 3*L/5, 4*L/5, 0, 0, m*g(0) * std::cos(std::atan(3/4));
    state.t     = 0;

    // Проинициализируем значение шага и количество иттераций
    double step = 0.01;
    unsigned iterations = 400;

    // Вызов метода численного решения
    std::vector<Vec> runge_kutta = implicitRK(state, step, iterations, rightPart, table);

    // Построение графиков полученных решений
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
    canvas.size(1500, 1000);
    canvas.save("runge_kutta.png");

    return 0;
}