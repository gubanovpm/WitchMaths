#include "../libs/shmemes.hh"

using namespace WitchMath;

int main(int argc, const char *argv[]) {
    matrix_t *matrix;
    get_args(argc, argv, &matrix);	
    matrix->compute();
    matrix->view_plot(0x42);
    while (true) {
        size_t plot_num;
        std::cout << std::endl << "What plot you want to see?" << std::endl;
        std::cin >> plot_num ;
        if (plot_num == (size_t)-1) break;
	    matrix->view_plot(plot_num);
    }
    delete matrix;
    return 0;
}