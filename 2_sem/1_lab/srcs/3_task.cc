#include "../libs/diff_eq.hh"

// Объявление констант
int M = 5;
double k = M_PI * M;

double beg_t = 0, end_t = 1;  // x ∈ [0, 1]
double y_0 = 0, y_1 = 0;

double L = end_t - beg_t;

// Функция правой части дифференциального уравнения y' = f(t, y) 
Vec rightPart(const Time& t, const Vec& s) noexcept {
    Vec result(s.size());
// Значение функций опредеяем здесь:
    double y = s(0), dy = s(1);
    // dy
    result(0) = dy ;
    // d2y
    result(1) = -std::pow(k, 2) * y ;

    return result;
}
//=======================================================================

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
    unsigned iterations = 100000;
    double step = (end_t - beg_t)/iterations;

    // Проинициализируем начальные значения
    double eps = 1e-6;
    Vec segment = Vec(2); segment << beg_t, end_t;
    Vec tg = Vec(2); tg << k, 3*k; // Пристрелочные углы наклона
    Vec state = Vec(2); state << y_0, y_1;

    // Вызов метода численного решения
    std::vector<Vec> res = shooting_method(segment, state, tg, rightPart, table, step, eps);

    // Построение графиков полученных решений
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.legend().atOutsideBottom().displayHorizontal().displayExpandWidthBy(2);

    sciplot::Vec x = sciplot::linspace(beg_t, end_t, iterations);
    sciplot::Vec y = sciplot::linspace(0, 0, iterations);
    for (unsigned i = 0; i < res.size(); ++i) { 
        y[i] = res[i](0); 
    }

    plot.drawCurve(x, y).label("y(x)").lineWidth(2);
    plot.grid().lineWidth(2).show();
    sciplot::Figure figure = {{ plot }};
    sciplot::Canvas canvas = {{ figure }};

    canvas.size(1000, 800);
    canvas.save("third_task.png");
    canvas.show();

    return 0;
}