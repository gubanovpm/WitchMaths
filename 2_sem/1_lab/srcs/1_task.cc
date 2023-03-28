#include "../libs/diff_eq.hh"
#include "../libs/tables.hh"

double g(const double t) {
    return 9.81 + 0.01 * std::cos(2 * M_PI * t);
    return 9.81;
}

double gt(const double t) {
    return -0.01 * std::sin(2*M_PI*t) * 2 * M_PI ;
    return 0;
}

double L = 4;
double m = 2;
/** Функция правой части дифференциального уравнения y' = f(t, y)
 * @param time время t
 * @param state состояние y
 * @return правая часть f(t, y)
 */ 
Vec rightPart(const Time& t, const Vec& s) {
    Vec result(s.size());
// Значение функций опредеяем здесь:
    double x = s(0), y = s(1), vx = s(2), vy = s(3), T = s(4);

    result(0) = vx;
    result(1) = vy;
    result(2) = -x/m/L*T;        // vx_t
    result(3) = -y/m/L*T + g(t); // vy_t
    result(4) = 2* m/L * (vx * result(2) + vy*result(3)) + m*gt(t)*y/L + m*g(t)*vy/L;

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
    table.column = std::array<double, 4> {   0, 1./2, 1./2,    1};
    table.string = std::array<double, 4> {1./6, 1./3, 1./3, 1./6};
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
    unsigned iterations = 1000;
    double beg_t = 0, end_t = 4;
    double step = (end_t - beg_t)/iterations;

    // Проинициализируем начальные значения
    double vx = 4;
    State state;
    state.state = Vec(5); state.state << 0, L, vx, 0, m/L*(vx*vx) + m*g(0) ;
    state.t     = beg_t;

    // Вызов метода численного решения
    const long unsigned stage = iterations; // стадийность для метода Адамса
    std::vector<Vec> runge_kutta = RungeKutta(state, step, stage, rightPart, table);
    std::vector<Vec> addams = explicitAdams(state, step, iterations, rightPart, coefs_3, runge_kutta);

    std::vector<Vec> runge_kutta_1 = RungeKutta(state, step, iterations, rightPart, table);


    // Построение графиков полученных решений
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.legend().atOutsideBottom().displayHorizontal().displayExpandWidthBy(2);

    sciplot::Vec t = sciplot::linspace(beg_t, end_t, iterations);
    sciplot::Vec x = sciplot::linspace(0, 0, iterations);
    sciplot::Vec y = sciplot::linspace(0, 0, iterations);
    sciplot::Vec x_1 = sciplot::linspace(0, 0, iterations);
    sciplot::Vec y_1 = sciplot::linspace(0, 0, iterations);
    for (unsigned i = 0; i < runge_kutta.size(); ++i) {
        x[i] = addams[i](0);
        y[i] = -addams[i](1);

        x_1[i] = runge_kutta_1[i](0);
        y_1[i] = -runge_kutta_1[i](1);
    }
    plot.drawCurve(x, y).label("z(x)").lineWidth(2);
    plot.drawCurve(x_1, y_1).label("z(x)").lineWidth(2);
    plot.grid().lineWidth(2).show();
    sciplot::Figure figure = {{ plot }};
    sciplot::Canvas canvas = {{ figure }};

    canvas.size(1000, 700);
    canvas.save("adams.png");
    canvas.show();

    return 0;
}