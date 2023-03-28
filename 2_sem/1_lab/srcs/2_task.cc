#include "../libs/diff_eq.hh"

// Объявление констант
// double m = 0.01215088;
double m = 0.012277471;
double M = 1 - m;
double f = 1;

// Объвление вспомогательных функций
double r_1(const double x, const double y) {
    return std::sqrt(std::pow(x + m, 2) + std::pow(y, 2));
}
double r_2(const double x, const double y) {
    return std::sqrt(std::pow(x - M, 2) + std::pow(y, 2));
}

/** Функция правой части дифференциального уравнения y' = f(t, y)
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */ 
Vec rightPart(const Time& t, const Vec& s) noexcept {
    Vec result(s.size());
// Значение функций опредеяем здесь:
    double x = s(0), y = s(1), vx = s(2), vy = s(3);
    // dx
    result(0) = vx ;
    // dy
    result(1) = vy ;
    // dvx
    result(2) = 2*vy + x - M*(x+m)/std::pow(r_1(x, y), 3) - m * (x - M)/std::pow(r_2(x, y), 3) - f * vx ;
    // dvy
    result(3) = -2*vx + y - M * y/std::pow(r_1(x, y), 3) - m * y/std::pow(r_2(x, y), 3) - f * vy;
    return result;
}

int main() {
    // Таблица Бутчера для явного метода Рунге-Кутты 2 порядка
    // ButcherTable<2> table;
    // table.column = std::array<double, 2> { 0, 1./2};
    // table.string = std::array<double, 2> { 0, 1.};
    // table.matrix = std::array<std::array<double, 2>, 2>  {
    //     std::array<double, 2> {    0,    0}, 
    //     std::array<double, 2> { 1./2,    0}, 
    // };

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
    unsigned iterations = 300000;
    double beg_t = 0, end_t = 16;
    double step = (end_t - beg_t)/iterations;

    // Проинициализируем начальные значения
    State state;
    double x_0  = 1.2, y_0 = 0, dx_0 = 0, dy_0 =  -1.05;
    // double x_0 = .994, y_0 = 0, dx_0 = 0, dy_0 = -2.031732629557337;
    state.state = Vec(4); state.state << x_0, y_0, dx_0, dy_0;
    state.t     = beg_t;

    // Вызов метода численного решения
    std::vector<Vec> runge_kutta = RungeKutta(state, step, iterations, rightPart, table);

    // Построение графиков полученных решений
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.legend().atOutsideBottom().displayHorizontal().displayExpandWidthBy(2);

    unsigned tt = 0, kk = 1;
    sciplot::Vec t = sciplot::linspace(beg_t, end_t, iterations/kk);
    sciplot::Vec x = sciplot::linspace(0, 0, iterations/kk);
    sciplot::Vec y = sciplot::linspace(0, 0, iterations/kk);
    for (unsigned i = 0; i < runge_kutta.size(); ++i) {
        if (i % kk == 0 ) {
            x[tt] = runge_kutta[i](0);
            y[tt] = -runge_kutta[i](1);
            ++tt;
        }
        // std::cout << x[i] << " ; " << y[i] << std::endl;
    }
    plot.drawCurve(x, y).label("z(x)").lineWidth(2);
    plot.grid().lineWidth(2).show();
    sciplot::Figure figure = {{ plot }};
    sciplot::Canvas canvas = {{ figure }};

    canvas.size(1000, 700);
    canvas.save("second_task.png");
    canvas.show();

    return 0;
}