#include <iostream>
#include <cmath>
#include <exception>
#include <functional>

#define L 4.
#define m 2.

#define eps 1e-4
#define CMP(__arg__) ((std::abs(__arg__) < eps))
//========================================================================
using ftype = std::function<double>(const double);
//========================================================================
class point_t {
public:
	double x;
	double y;

	point_t(const double x, const double y): x(x), y(y) {}

	void dump(std::ostream &out) const {
		out << "(" << x << ", " << y << ")" ; 
	}
};

std::ostream &operator<<(std::ostream &out, point_t &p);

double get_ocoord(const point_t &p);
//========================================================================
int main() {
	point_t p1(1, 0);
	std::cout << p1 << std::endl;
	std::cout << get_ocoord(p1) << std::endl;
	return 0;
}
//========================================================================

std::ostream &operator<<(std::ostream &out, point_t &p) { 
	p.dump(out);
	return out;
}

double get_ocoord(const point_t &p) {
	if (CMP(p.x)) { return L*L - p.y * p.y; }
	if (CMP(p.y)) { return L*L - p.x * p.x; }
	throw std::invalid_argument("Wrong coordinates\n");
}
//========================================================================
