#include "../libs/schemes.hh"
using namespace WitchMath;

double border_func(const double x) {
    //return 19 * std::exp(-10 * std::pow((x - 50), 2)) + 1;
	return ((x >= 0.25 && x <= 0.75) ? 0.5 : 0);
}

int main(int argc, const char *argv[]) {
	matrix_t *matrix;
	get_args(argc, argv, &matrix, border_func);	
	matrix->compute();
	matrix->save_img("../imgs/test/test", "lux.gif");
	delete matrix;
	return 0;
}