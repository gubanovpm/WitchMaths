#include "../libs/schemes.hh"

using namespace WitchMath;

void WitchMath::get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u) noexcept {
	if (argc != 8) {
		std::cout << "Usage is: " << argv[0] << " <T> <X> <k> <m> <alpha> <beta>" << std::endl;
		std::cout << "<T> -- max X value (min X = 0)" << std::endl;
		std::cout << "<X> -- max T value (min T = 0)" << std::endl;
		std::cout << "<k> -- count of steps on X coord" << std::endl;
		std::cout << "<m> -- count of steps on T coord" << std::endl;
        std::cout << "<alpha> -- alpha coefficient" << std::endl;
        std::cout << "<beta> -- beta coefficient" << std::endl;
        std::cout << "<type> -- scheme type" << std::endl;
		exit(-1);
	}
	*matrix = new matrix_t {
		std::atof(argv[1]), 
		std::atof(argv[2]), 
		std::atoi(argv[3]), 
		std::atoi(argv[4]),
		u,
        std::atof(argv[5]),
        std::atof(argv[6]), 
        std::atoi(argv[7])
    };
}

void matrix_t::compute() {	
	double h = x;
	for (size_t j = 0; j < m; ++j) { data[j] = u(j*h); }


}

double matrix_t::scheme_1(const size_t i, const size_t j, const double c) {
    size_t indn = ((j   == 0) ? j : j - 1);
    size_t indx = ((j+1 == m) ? m : j + 1);
    double coef = c*t/2/x/x;
    return coef*(data[(i+1)*m + indx]-2*data[])
    ;
}

void matrix_t::save_img(const std::string path, const std::string gif_name) {
	// std::system("rm ../imgs/lux_scheme/*.png");
	// for (size_t i = 0; i < n; ++i) {
	// 	auto f = matplot::figure(true);
	// 	auto J = matplot::linspace(0, X, m);
	// 	std::vector<double> U;
	// 	for (size_t j = 0; j < m; ++j) U.push_back(data[i*m + j]);
	// 	matplot::plot(J, U, "--xr");
	// 	matplot::xlabel("x");
	// 	matplot::ylabel("U in " +  std::to_string(t*i));
	// 	matplot::save(path + "_" + std::to_string(t*i)+".png");
	// }
	// std::system("convert -delay 1 -loop 0 ../imgs/temp/*.png lux.gif");
	// std::system("rm ../imgs/lux_scheme/*.png");
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
