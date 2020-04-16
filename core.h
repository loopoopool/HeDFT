#ifndef CORE_H
#define CORE_H
#include <cmath>
#include <functional>
#include <Eigen/Dense>

using namespace Eigen;

const double EPS = 1e-15;

double simpson(ArrayXd u, double step){
	double a = 0.;
	for (int i=0; i < (u.size()-2)/2.; ++i)
		a += step*(u.coeff(2*i) + 4.*u.coeff(2*i+1) + u.coeff(2*i+2))/3.;
	return a;
}

template <typename X, typename ...T>
using Func = std::function<double(X&, double&, T&...)>;

template <typename X, typename ...T> ArrayXd VerletREV(double y0, double y1, 
		Eigen::Array<X, Dynamic, 1> &t, double step, Func<X, T...> &f, T& ...par){
	int n = t.size();
	ArrayXd y(n);
	y(n-1) = y0;
	y(n-2) = y1;
	for (int i=n-2; i>0; --i)
		y(i-1) = 2.*y(i) - y(i+1) + pow(step, 2)*f(t(i), y(i), par...);
	return y;
}

template <typename X, typename ...T> ArrayXd Verlet(double y0, double y1, 
		Eigen::Array<X, Dynamic, 1> &t, double step, Func<X, T...> &f, T& ...par){
	int n = t.size();
	ArrayXd y(n);
	y[0] = y0;
	y[1] = y1;
	for (int i=1; i<n-1; ++i)
		y(i+1) = 2.*y(i) - y(i-1) + pow(step, 2)*f(t(i), y(i), par...);
	return y;
}

#endif
