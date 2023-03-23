#include "../libs/diff_eq.hh"

// Объявление констант
const double mu    = 1./82.45;
const double gamma = 1 - mu;
const double f     = 0;

// Объвление вспомогательных функций
double r_1(const double x, const double y) {
    return std::sqrt(std::pow(x+mu, 2) + std::pow(y, 2));
}
double r_2(const double x, const double y) {
    return std::sqrt(std::pow(x-gamma, 2) + std::pow(y, 2));
}

/** Функция правой части дифференциального уравнения y' = f(t, y)
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */ 
Vec rightPart(const Time& t, const Vec& s) noexcept {
    Vec result(s.size());
// Значение функций опредеяем здесь:
    
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

    // Проинициализируем значение шага и количество иттераций
    unsigned iterations = 8000;
    double beg_t = 0, end_t = 8;
    double step = (end_t - beg_t)/iterations;

    // Проинициализируем начальные значения
    State state;
    double x_0  = 1.2, y_0 = -1.05, dx_0 = 0;
    double dy_0 = (-2*dx_0 + y_0 - gamma * y_0/std::pow(r_1(x_0, y_0), 3) - mu * y_0/std::pow(r_2(x_0, y_0), 3) ) / (1+f);
    state.state = Vec(4); state.state << x_0, y_0, dx_0,  dy_0;
    state.t     = beg_t;

    // Вызов метода численного решения
    std::vector<Vec> runge_kutta = implicitRK(state, step, iterations, rightPart, table);

    // Построение графиков полученных решений
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.legend().atOutsideBottom().displayHorizontal().displayExpandWidthBy(2);

    sciplot::Vec t = sciplot::linspace(beg_t, end_t, iterations);
    sciplot::Vec x = sciplot::linspace(0, 0, iterations);
    sciplot::Vec y = sciplot::linspace(0, 0, iterations);
    for (unsigned i = 0; i < runge_kutta.size(); ++i) {
        x[i] = runge_kutta[i](0);
        y[i] = -runge_kutta[i](1);
        // std::cout << x[i] << " ; " << y[i] << std::endl;
    }
    plot.drawCurve(x, y).label("z(x)").lineWidth(2);
    plot.grid().lineWidth(2).show();
    sciplot::Figure figure = {{ plot }};
    sciplot::Canvas canvas = {{ figure }};

    canvas.size(1000, 700);
    canvas.save("runge_kutta.png");
    canvas.show();

    return 0;
}