#include "schemes.hh"

using namespace WitchMath;

void WitchMath::get_args(const int argc, const char *argv[], matrix_t **matrix, const utype u) noexcept {
	if (argc != 7) {
		std::cout << "Usage is: " << argv[0] << " <X> <T> <m> <k> <type> <is_div>" << std::endl;
		std::cout << "<X> -- max X value (min X = 0)" << std::endl;
		std::cout << "<T> -- max T value (min T = 0)" << std::endl;
		std::cout << "<m> -- count of steps on X coord" << std::endl;
		std::cout << "<k> -- count of steps on T coord" << std::endl;
		std::cout << "<type> -- using scheme type:" << std::endl;
		std::cout << "\t0 - Lux scheme" << std::endl;
		std::cout << "\t1 - Lux-Wendorf scheme" << std::endl;
		std::cout << "\t2 - Right corner scheme" << std::endl;
		std::cout << "\t3 - Left corner scheme" << std::endl;
		std::cout << "<is_div> - is div form" << std::endl;
		exit(-1);
	}
	*matrix = new matrix_t {
		std::atof(argv[1]), 
		std::atof(argv[2]), 
		std::atoi(argv[3]), 
		std::atoi(argv[4]),
		std::atoi(argv[5]), 
		u,
		std::atoi(argv[6])};	
}

void matrix_t::compute() {	
	double h = x;
	for (size_t j = 0; j < m; ++j) {
		data[j] = u(j*h);
		// std::cout << data[j] << " " ;
	}
	std::cout << std::endl;
	switch (type) {
		case 0: { //Lux scheme
			if (!is_div) { 
				std::cout << "Lux nodiv" << std::endl;
				for (size_t i = 1; i < n; ++i) 
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = Lux_scheme_nodiv(i-1, j, data[(i-1)*m + j]);
					}
			} else {
				std::cout << "Lux div" << std::endl;
				for (size_t i = 1; i < n; ++i) 
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = Lux_scheme_div(i-1, j, 1);
					}
			}
			break;
		}
		case 1: {
			if (!is_div) { 
				std::cout << "Lux-Wendorf nodiv" << std::endl;
				for (size_t i = 1; i < n; ++i) 
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = Lux_Wendorf_scheme_nodiv(i-1, j, data[(i-1)*m + j]);
					}
			} else { 
				std::cout << "Lux-Wendorf div" << std::endl;
				for (size_t i = 1; i < n; ++i)
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = Lux_Wendorf_scheme_div(i-1, j, 1);
					}
			}
			break;
		}
		case 2: {
			if (!is_div) {
				std::cout << "Right corner nodiv" << std::endl;
				for (size_t i = 1; i < n; ++i)
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = right_corner_scheme_nodiv(i-1, j, data[(i-1)*m + j]);
					}
			} else { 
				std::cout << "Right corner div" << std::endl;
				for (size_t i = 1; i < n; ++i)
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = right_corner_scheme_div(i-1, j, 1);
					}
			}
			break;
		}
		case 3: {
			if (!is_div) {
				std::cout << "Left corner nodiv" << std::endl;
				for (size_t i = 1; i < n; ++i) {
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = left_corner_scheme_nodiv(i-1, j, data[(i-1)*m + j]);
						// std::cout << data[i*m + j] << " " ;
					}
					// std::cout << std::endl;
				}

			} else {
				std::cout << "Left corner div" << std::endl;
				for (size_t i = 1; i < n; ++i)
					for (size_t j = 0; j < m; ++j) {
						data[i*m + j] = left_corner_scheme_div(i-1, j, 1);
					}
			}
			break;
		}
		default: {
			std::cout << "Unknown scheme type" << std::endl;
			exit(-1);
		}
	}
}

double matrix_t::Lux_scheme_nodiv(const size_t i, const size_t j, const double c) const {
	double h = x;
	if (j   == 0) return (data[i*m+j+1]*(h-c*t) + data[i*m+  j]*(h+c*t)) / (2*h);
	if (j+1 == m) return (data[i*m+  j]*(h-c*t) + data[i*m+j-1]*(h+c*t)) / (2*h);
	return               (data[i*m+j+1]*(h-c*t) + data[i*m+j-1]*(h+c*t)) / (2*h);
}

double matrix_t::Lux_scheme_div(const size_t i, const size_t j, const double c) const {
	double h = x;
	if (j   == 0) return (data[i*m+j+1]*(h-data[i*m+j+1]*c*t/2) + data[i*m+  j]*(h+data[i*m+  j]*c*t/2)) / (2*h);
	if (j+1 == m) return (data[i*m+  j]*(h-data[i*m+  j]*c*t/2) + data[i*m+j-1]*(h+data[i*m+j-1]*c*t/2)) / (2*h);
	return               (data[i*m+j+1]*(h-data[i*m+j+1]*c*t/2) + data[i*m+j-1]*(h+data[i*m+j-1]*c*t/2)) / (2*h);
}

double matrix_t::Lux_Wendorf_scheme_nodiv(const size_t i, const size_t j, const double c) const {
	double h = x;
	double c1 = c*t/2/h;
	double c2 = c1*t/h;
	double prev = ((j   == 0) ? data[i*m+j] : data[i*m+j-1]);
	double next = ((j+1 == m) ? data[i*m+j] : data[i*m+j+1]);
	return data[i*m+j]*(1-2*c2) + next*(c2-c1) + prev*(c2+c1); 
}
double matrix_t::Lux_Wendorf_scheme_div(const size_t i, const size_t j, const double c) const {
	double h = x;
	double c1 = c*t/2/h;
	double c2 = c1*t/h;
	double prev = ((j   == 0) ? data[i*m+j] : data[i*m+j-1]);
	double next = ((j+1 == m) ? data[i*m+j] : data[i*m+j+1]);
	return data[i*m+j]*(1-data[i*m+j]*c2) + next*next*(c2-c1)/2 + prev*prev*(c2+c1)/2; 
}

double matrix_t::right_corner_scheme_nodiv(const size_t i, const size_t j, const double c) const{
	double h = x;
	double k1 = c*t/h;
	double prev = ((j == 0) ? data[i*m+j] : data[i*m+j-1] );
	return data[i*m+j]*(1-k1)+k1*prev;
}
double matrix_t::right_corner_scheme_div(const size_t i, const size_t j, const double c) const {
	double h = x;
	double k1 = c*t/h;
	double prev = ((j == 0) ? data[i*m+j] : data[i*m+j-1] );
	return data[i*m+j]*(1-k1*data[i*m+j]/2)+k1*prev*prev/2;
}

double matrix_t::left_corner_scheme_nodiv(const size_t i, const size_t j, const double c) const {
	double h = x;
	double k1 = c*t/h;
	double next = ((j +1 == m) ? data[i*m+j] : data[i*m+j+1] );
	return data[i*m+j]*(1.+c*k1)-k1*next;
}
double matrix_t::left_corner_scheme_div(const size_t i, const size_t j, const double c) const{
	double h = x;
	double k1 = c*t/h;
	double next = ((j +1 == m) ? data[i*m+j] : data[i*m+j+1] );
	return data[i*m+j]*(1.+c*k1*data[i*m+j]/2)-k1*next*next/2;
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
	// std::system("convert -delay 1 -loop 0 ../imgs/lux_scheme/*.png lux.gif");
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
