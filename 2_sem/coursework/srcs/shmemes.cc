#include "../libs/shmemes.hh"

using namespace WitchMath;
//===============================================================================================
//===============================================================================================
void WitchMath::get_args(const int argc, const char *argv[], matrix_t **matrix) noexcept {
	if (argc != 5) {
		std::cout << "Usage is: " << argv[0] << " <Tmax> <Rmax> <k> <m>" << std::endl;
		std::cout << "<Tmax> -- max Tmax value (Tmin = 0 by default)" << std::endl;
		std::cout << "<Rmax> -- max Rmax value (Rmin = 0 by default)" << std::endl;
		std::cout << "<k>    -- count of steps on X coord" << std::endl;
		std::cout << "<m>    -- count of steps on T coord" << std::endl;
		exit(-1);
	}
	*matrix = new matrix_t {
		std::atof(argv[1]), 
		std::atof(argv[2]), 
		std::atoi(argv[3]), 
		std::atoi(argv[4]),
    };
}
//===============================================================================================
//===============================================================================================
void matrix_t::scheme_u (const size_t i, const size_t j, const double c) noexcept {
    double k = c * tau / h;
    size_t jmn = ((j  == 0) ? j : j-1);
    size_t jmx = ((j+1== m) ? j : j+1);

    u[i+1][j] = k/2 * (
        (k - 1) * P[i][jmx] -
        2 * P[i][j] +
        (k + 1) * P[i][jmn]
    ) + u[i][m];
}
void matrix_t::scheme_u (const size_t i, const double c) noexcept {
    double k = c * tau / h;
    u[i+1][0] = k/2 * (k - 1) * P[i][1]+ u[i][0];
}
//===============================================================================================
//===============================================================================================
void matrix_t::scheme_ro(const size_t i, const size_t j, const double c) noexcept {
    double k = c * tau / h;
    size_t jmn = ((j  == 0) ? j : j-1);
    size_t jmx = ((j+1== m) ? j : j+1);
    ro[i+1][j] = k/2 * (
        (k - 1) * std::pow(jmx*h, 2) * ro[i][jmx] * u[i][jmx] - 
        2 * std::pow(j*h, 2) * ro[i][j] * u[i][j] +
        (k + 1) * std::pow(jmn*h, 2) * ro[i][jmn] * u[i][jmn]
    ) + ro[i][j];
}
void matrix_t::scheme_ro(const size_t i, const double c) noexcept {
    double k = c * tau / h;
    ro[i+1][0] = k/2 * (k - 1) * std::pow(h, 2) * ro[i][1] * u[i][1] + ro[i][0];
}
//===============================================================================================
//===============================================================================================
void matrix_t::scheme_EP(const size_t i, const size_t j, const double c) noexcept {
    double k = c * tau / h, result = 0;
    size_t jmn = ((j  == 0) ? j : j-1);
    size_t jmx = ((j+1== m) ? j : j+1);

    E[i+1][j] = k/2 * (
        (k - 1) * std::pow(jmx*h, 2) * u[i][jmx] -
        2 * std::pow(j*h, 2) * u[i][j] +
        (k + 1) * std::pow(jmn*h, 2) * u[i][jmn]
    ) + E[i][j];
    P[i+1][j] = E[i+1][j] * ro[i+1][j] * (g - 1);
}
void matrix_t::scheme_EP(const size_t i, const double c) noexcept {
    double k = c * tau / h, result = 0;
    E[i+1][0] = k/2 * (k - 1) * std::pow(h, 2) * u[i][1] + E[i][0];
    P[i+1][0] = E[i+1][0] * ro[i+1][0] * (g - 1);
}
//===============================================================================================
//===============================================================================================
void matrix_t::compute() {	
    for (size_t i = 0, end = n-1; i < end; ++i) {
        scheme_ro(i, 0, 1/(std::pow(h, 2)));
        scheme_u (i, 0, 1/ro[i][0]);
        scheme_EP(i, 0, P[i][0]/(std::pow(h, 2)*ro[i][0]));

        for (size_t j = 1; j < m; ++j) {
            scheme_ro(i, j, 1/(std::pow(j*h, 2)));
            if (ro[i][j] == 0) { break; }

            scheme_u (i, j, 1/ro[i][j]);
            scheme_EP(i, j, P[i][j]/(std::pow(j*h, 2)*ro[i][j]));
        }
    }
}
//===============================================================================================
//===============================================================================================
void matrix_t::view_plot(const size_t plot_num) noexcept {
    std::vector<std::vector<double>> *data = nullptr;
    switch (plot_num) {
        case 0: { data = &ro; break;}
        case 1: { data =  &u; break;}
        case 2: { data =  &E; break;}
        case 3: { data =  &P; break;}
        default: {
            std::cout << "Plot types is: " << std::endl;
            std::cout << "\t 0 -- ro(t, r)" << std::endl;
            std::cout << "\t 1 -- u(t, r)" << std::endl;
            std::cout << "\t 2 -- E(t, r)" << std::endl;
            std::cout << "\t 3 -- P(t, r)" << std::endl;
            std::cout << "\t-1 -- if you want to stop programm" << std::endl;
            return;
        }
    }

	auto f = matplot::figure(true);
	auto [R, Y] = matplot::meshgrid(
                    matplot::linspace(0, T, n), 
                    matplot::linspace(0, X, m)
                );
	auto [I, J] = matplot::meshgrid(
                    matplot::linspace(0, n-1, n), 
                    matplot::linspace(0, m-1, m)
                );
	auto Z = matplot::transform(I, J, [data](size_t i, size_t j) {  
		return (*data)[i][j];
	});

	matplot::contour(Y, R, Z);
    matplot::surf(Y, R, Z);
	matplot::xlabel("r");
	matplot::ylabel("t");
	//matplot::save(path);
	matplot::show();
}
