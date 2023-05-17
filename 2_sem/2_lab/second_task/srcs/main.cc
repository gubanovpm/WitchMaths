#include "../libs/schemes.hh"
using namespace WitchMath;

#define XMX 20

double border_func(const double x) {
	return ((x >= 0.25 * XMX && x <= 0.75 * XMX) ? 0.5 : 0);
}
double g(const double t) {
	return t;
}

int main(int argc, const char *argv[]) {
	matrix_t *matrix;
	get_args(argc, argv, &matrix, border_func, g);	
	matrix->compute();
	matrix->save_img();
	delete matrix;
	return 0;
}