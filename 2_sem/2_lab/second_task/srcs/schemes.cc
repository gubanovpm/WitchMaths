#include "../libs/schemes.hh"

using namespace WitchMath;

void WitchMath::get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u, const utype g) noexcept {
	if (argc != 7) {
		std::cout << "Usage is: " << argv[0] << " <T> <X> <k> <m> <alpha> <beta>" << std::endl;
		std::cout << "<T> -- max X value (min X = 0)" << std::endl;
		std::cout << "<X> -- max T value (min T = 0)" << std::endl;
		std::cout << "<k> -- count of steps on X coord" << std::endl;
		std::cout << "<m> -- count of steps on T coord" << std::endl;
        std::cout << "<alpha> -- alpha coefficient" << std::endl;
        std::cout << "<beta> -- beta coefficient" << std::endl;
		exit(-1);
	}
	*matrix = new matrix_t {
		std::atof(argv[1]), 
		std::atof(argv[2]), 
		std::atoi(argv[3]), 
		std::atoi(argv[4]),
		u,
		g,
        std::atof(argv[5]),
        std::atof(argv[6])
    };
}

void matrix_t::compute() {	
	double h = x;
	for (size_t j = 0; j < m; ++j) { data[0][j] = u(j*h); }
	for (size_t i = 0; i < n; ++i) { data[i][0] = g(i*t); }

	for (size_t i = 1; i < n; ++i) {
		for (size_t j = 1; j < m; ++j) {
			data[i][j] = scheme_1(i - 1, j, 1);
		}
	}
}

double matrix_t::scheme_1(const size_t i, const size_t j, const double c) {
	// std::cout << i << " ; " << j << std::endl;
    size_t jmn = ((j   == 0) ? j : j - 1);
    size_t jmx = ((j+1 == m) ? m : j + 1);
	size_t imn = ((i   == 0) ? i : i - 1);
    double amin = amn(i, jmn);
	double amax = amn(i, j);
    return (
		data[i+1][j] * (std::pow(x, 2) + amax + amin) +
		data[imn][jmn] * amin -
		data[i][j] * std::pow(x, 2) / t - 
		std::pow(x, 2) * Q(i*t)
	) / amn(i, j);
    ;
}

void matrix_t::save_img() noexcept{
	auto f = matplot::figure(true);
	auto [R, Y] = matplot::meshgrid(matplot::linspace(0, T, n), matplot::linspace(0, X, m));
	auto [I, J] = matplot::meshgrid(matplot::linspace(0, n-1, n), matplot::linspace(0, m-1, m));
	auto Z = matplot::transform(I, J, [this](size_t i, size_t j) {  
		return data[i][j];
	});

	matplot::surf(Y, R, Z);
	matplot::xlabel("x");
	matplot::ylabel("t");
	//matplot::save(path);
	matplot::show();
}
