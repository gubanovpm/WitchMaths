#ifndef __newton_hh__
#define __newton_hh__

#include <iostream> 
#include <vector>

namespace witch_math {

struct newton_interpolation {
private:
	std::vector<double> _coefs = {};
	std::vector<double> _x = {};
public:
	newton_interpolation(const std::vector<double> &, const std::vector<double> &);
	void dump (std::ostream &) const;
	double evaluate(double) const;
	void evaluate(const std::vector<double> &, std::vector<double> &) const;

};


}

std::ostream &operator<< (std::ostream &, witch_math::newton_interpolation &);
std::ostream &operator<< (std::ostream &, const witch_math::newton_interpolation &);

#endif

