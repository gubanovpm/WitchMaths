#include "schemes.hh"

using namespace WitchMath;

void WitchMath::get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u) noexcept {
	if (argc != 6) {
		std::cout << "Usage is: " << argv[0] << " <X> <T> <m> <k> <type>" << std::endl;
		std::cout << "<X> -- max X value (min X = 0)" << std::endl;
		std::cout << "<T> -- max T value (min T = 0)" << std::endl;
		std::cout << "<m> -- count of steps on X coord" << std::endl;
		std::cout << "<k> -- count of steps on T coord" << std::endl;
		std::cout << "<type> -- using scheme type:" << std::endl;
		std::cout << "\t0 - Lux scheme" << std::endl;
		std::cout << "\t1 - Lux-Wendorf scheme" << std::endl;
		std::cout << "\t2 - Right corner scheme" << std::endl;
		std::cout << "\t3 - Left corner scheme" << std::endl;
		exit(-1);
	}
	*matrix = new matrix_t {std::atof(argv[1]), std::atof(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]), u};	
}

void matrix_t::compute() {	
	double h = x;
	for (size_t j = 0; j < m; ++j) {	
		data[j] = u(j*h);
	}
	switch (type) {
		case 0: { //Lux scheme
			printf("starting computing Lux scheme, with t = %lg, h = %lg\n", t, h);
			for (size_t i = 1; i < n; ++i)
				for (size_t j = 0; j < m; ++j) {
					Lux_scheme(i-1, j, data[(i-1)*m + j]);
				}
		}
	}
}

double matrix_t::Lux_scheme(const size_t i, const size_t j, const double c) const {
	double h = x;
	if (j   == 0) return (data[i*m+j+1]*(h-c*t) + data[i*m+  j]*(h+c*t)) / (2*h);
	if (j+1 == m) return (data[i*m+  j]*(h-c*t) + data[i*m+j-1]*(h+c*t)) / (2*h);
	return               (data[i*m+j+1]*(h-c*t) + data[i*m+j-1]*(h+c*t)) / (2*h);
}

void matrix_t::save_img(const std::string path) {
	auto f = matplot::figure(true);
	auto [R, Y] = matplot::meshgrid(matplot::linspace(0, T, n), matplot::linspace(0, X, m));
	auto [I, J] = matplot::meshgrid(matplot::linspace(0, n-1, n), matplot::linspace(0, m-1, m));
	auto Z = matplot::transform(I, J, [this](size_t i, size_t j) {  
		return data[i*m + j];
	});

	matplot::surf(Y, R, Z);
	matplot::xlabel("X");
	matplot::ylabel("T");
	//matplot::save(path);
	matplot::show();
}
