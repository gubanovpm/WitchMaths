#include "../libs/integral.hh"

int main() {
	try {
		auto f = [](const double x){ return std::sin(100*x) * std::exp(-std::pow(x, 2))*std::cos(2*x); };
		auto g = [](const double x){ return -x/std::tan(x); };

		double f_right_answer = 0.010006097860331;
		double g_right_answer = -1.08879;

		std::vector<witch_math::IIntegral *> integral_f = {};
		std::vector<std::string> names = {"Rectangle: ", "Trapezoid: ", "Simpsons: ", "3/8: "};
		integral_f.push_back(new witch_math::rectangle_integral (f));
		integral_f.push_back(new witch_math::trapezoid_integral (f));
		integral_f.push_back(new witch_math::simpsons_integral (f));
		integral_f.push_back(new witch_math::three_eighths_integral (f));

		std::vector <double> x;
		size_t count = 500;
		for (size_t i = 0; i < count; ++i) {
			x.push_back(double(i) * 3. / count);
		}

		for (size_t i = 0; i < integral_f.size(); ++i) {
			double t_r = integral_f[i]->value(x);
			std::cout << names[i] << t_r << " with error = " << std::abs(f_right_answer - t_r) << std::endl;
		}
	
		std::cout << std::endl;

		x = {};
		double t = 0.01;
		for (auto el: integral_f) el->set_function(g);
		for (size_t i = 0; i < count; ++i) {
			x.push_back(t + double(i) * (M_PI - t) / 2 / count);
		}

		for (size_t i = 0; i < integral_f.size(); ++i) {
			double t_r = std::sin(t) + integral_f[i]->value(x);
			std::cout << names[i] << t_r << " with error = " << std::abs(g_right_answer - t_r) << std::endl;
		}


		for (size_t i = 0, end = integral_f.size(); i < end; ++i) {
			delete(integral_f[i]);
		}
		
	}
	catch(std::exception &e) {
		std::cout << e.what() << std::endl;
	}

	return 0;
}
