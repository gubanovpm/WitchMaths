#include "../libs/diff_eq.hh"
#include "../libs/tables.hh"

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
    // Проинициализируем значение шага и количество иттераций
    unsigned iterations = 100000;
    double step = (end_t - beg_t)/iterations;

    // Проинициализируем начальные значения
    double eps = 1e-6;
    Vec segment = Vec(2); segment << beg_t, end_t;
    Vec tg = Vec(2); tg << k, 3*k; // Пристрелочные углы наклона
    Vec state = Vec(2); state << y_0, y_1;

    // Вызов метода численного решения
    std::vector<Vec> res = shooting_method(segment, state, tg, rightPart, ex_table_4, step, eps);

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