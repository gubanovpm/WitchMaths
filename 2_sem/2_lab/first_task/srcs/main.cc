#include "schemes.hh"
using namespace WitchMath;

matrix_t *matrix;
double border_func(const double x) {
	return ((x >= 0.25 && x <= 0.75) ? 1 : 0);
}

int main() {
	//SchemeT t;
	size_t n, m, type;
	int res;

	std::cin >> n >> m >> type ;
	matrix = new matrix_t {n , m, border_func};

	pthread_t thid[THREAD_COUNT];
	thread_arg_t args[THREAD_COUNT];
	for (size_t i = 0; i < THREAD_COUNT; ++i) {
		args[i] = {i, THREAD_COUNT, type};
	}

	for (size_t i = 0; i < THREAD_COUNT; ++i) {
		if (res = pthread_create(thid + i, (pthread_attr_t *)NULL, compute, 0)) {
			std::cout << "Error on thread create" << std::endl;
			exit(-1);
		}
	}
	#ifdef HAVE_THR_SETCONCURRENCY_PROTO
		thr_setconcurrency(THREAD_COUNT);
	#endif

	for (size_t i = 0; i < THREAD_COUNT; ++i) {
		pthread_join(thid[i], (void **)NULL);
	}

	delete matrix;
	return 0;
}
