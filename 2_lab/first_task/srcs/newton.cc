#include "../libs/newton.hh"

witch_math::newton_interpolation::newton_interpolation(const std::vector<double> &x, const std::vector<double> &y) {
	if (x.size() != y.size()) {
		throw std::logic_error("Differnet size of Argument and Value vectors!");
	}

	std::vector<std::vector<double>> temporary = {{}};
	_coefs.push_back(y[0]);

	for (size_t i = 0; i < y.size(); ++i) {
		temporary[0].push_back(y[i]);
		_x.push_back(x[i]);
	}

	for (size_t layer = 1; layer < x.size() ; ++layer) {
		temporary.push_back({});
		for (size_t i = 1; i <= x.size() - layer; ++i) {
			temporary[layer].push_back((temporary[layer - 1][i] - temporary[layer - 1][i - 1]) / (x[layer] - x[0]));
		}
		_coefs.push_back(temporary[layer][0]);
	}

}

void witch_math::newton_interpolation::dump (std::ostream &out) const {
	for (size_t i = 0; i < _coefs.size() ; ++i) {
		if (_coefs[i] > 0)	out << " + " << _coefs[i];
		else								out << " - " << -1 * _coefs[i];

		for (size_t j = 0; j < i; ++j) {
			out << "(x" ;
			if (_x[j] > 0)	out << " - " << _x[j] ;
			else						out << " + " << -1 * _x[j] ;
			out << ")" ;
		}
	}
}

std::ostream &operator<< (std::ostream &out, witch_math::newton_interpolation &P) {
	P.dump(out);
	return out;
}

std::ostream &operator<< (std::ostream &out, const witch_math::newton_interpolation &P) {
	P.dump(out);
	return out;
}


double witch_math::newton_interpolation::evaluate(const double x) const {
	double result = _coefs[0];
	for (size_t i = 1; i < _coefs.size(); ++i) {
		double coef = 1;
		for (size_t j = 0; j < i; ++j) {
			coef *= (x - _x[j]);
		}
		result += coef * _coefs[i];
	}
	return result;
}

void witch_math::newton_interpolation::evaluate(const std::vector<double> &x, std::vector<double> &y) const {
	y = {};
	for (size_t i = 0; i < x.size(); ++i)
		y.push_back(evaluate(x[i]));
}
