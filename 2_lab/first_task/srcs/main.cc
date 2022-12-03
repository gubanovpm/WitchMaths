#include "../libs/newton.hh"

#include <matplot/matplot.h>

using namespace matplot;

int main() {
	
	const std::vector<double> t1 = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
	const std::vector<double> x = {1, 0.8, 0.5, 0.307, 0.2, 0.137, 0.1, 0.075, 0.06, 0.047, 0.039};

	const std::vector<double> t2 = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
	const std::vector<double> y = {0.02, 0.079, 0.175, 0.303, 0.459, 0.638, 0.831, 1.03, 1.23, 1.42};

	try {
		witch_math::newton_interpolation P(t1, x);
		witch_math::newton_interpolation M(t2, y);

		std::cout << "Pn = " << P << std::endl;
		std::cout << "Mn = " << M << std::endl;

		std::vector<double> t_1 = {}, f_1 = {}, t_2 = {}, f_2 = {};
		size_t points_count = 1000;
		for (size_t i = 0; i < points_count; ++i) {
			t_2.push_back(-0.8 + double(i) * 1.8/points_count); 
			t_1.push_back(double(i) * 1./points_count);
		}

		P.evaluate(t_1, f_1);
		M.evaluate(t_2, f_2);

		grid(true);
		plot(t_1, f_1)->line_width(1).color("red");
		show();

		grid(true);
		plot(t_2, f_2)->line_width(1).color("blue");
		show();

		M.evaluate(t_1, f_2);
		plot(f_1, f_2)->line_width(1).color("green");
		show();

	}
	catch (std::exception &e) {
		std::cout << e.what() << std::endl;
	}

	return 0;

}
