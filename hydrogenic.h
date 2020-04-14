#ifndef HYDROGENIC_H
#define HYDROGENIC_H

#include <iostream>
#include <cmath>
#include <functional>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double EPS = 1e-15;

double simpson(ArrayXd u, double step){
	double a = 0.;
	for (int i=0; i < (u.size()-2)/2.; ++i)
		a += step*(u.coeff(2*i) + 4.*u.coeff(2*i+1) + u.coeff(2*i+2))/3.;
	return a;
}

enum class Mode{ FORWARD, BACKWARD };

enum GridType{ HOMOGENEOUS, THIJSSEN };

template <typename X, typename ...T>
using Func = function<double(X&, double&, T&...)>;

template <typename X, typename ...T> ArrayXd Verlet(double y0, double y1, 
		Eigen::Array<X, Dynamic, 1> &t, double step, Mode m, Func<X, T...> &f, T& ...par){
	int n = t.size();
	ArrayXd y(n);
	Array<X, Dynamic, 1> tt(n);
	switch(m){
		case Mode::FORWARD:
			y(0) = y0;
			y(1) = y1;
			tt = t;
			break;
		case Mode::BACKWARD:
			y(n-1) = y0;
			y(n-2) = y1;
			y = y.reverse();
			tt = t.reverse();
			break;
		default: throw;
	}
	for (int i=1; i<n-1; ++i)
		y(i+1) = 2.*y(i) - y(i-1) + pow(step, 2)*f(tt(i), y(i), par...);
	if (m==Mode::BACKWARD) return y.reverse();
	else return y;
}

/*************************************************
 * CLASSES DECLARATIONS
 *************************************************/

struct Grid{
	Grid(unsigned long, double);
	Grid(const Grid &);
	virtual ~Grid() = default;
	virtual ArrayXd jacobian() const = 0;
	virtual ArrayXd wfFactor() const = 0;
	unsigned long nmax;
	double rmax, h;
	ArrayXd r;
};

struct HomogeneousGrid : Grid{
	HomogeneousGrid(unsigned long, double);
	HomogeneousGrid(const HomogeneousGrid&);
	~HomogeneousGrid() = default;
	ArrayXd jacobian() const;
	ArrayXd wfFactor() const;
};

struct ThijssenGrid : Grid{
	ThijssenGrid(unsigned long, double, double);
	ThijssenGrid(const ThijssenGrid&);
	~ThijssenGrid() = default;
	ArrayXd jacobian() const;
	ArrayXd wfFactor() const;
	double delta, rp;
};


template <typename GRID, typename ...T> class Verlet_HZorb{
	public:
		Verlet_HZorb(int, int, int);
		Verlet_HZorb(const Verlet_HZorb &);
		~Verlet_HZorb();
		ArrayXd u() const;
		ArrayXd r() const;
		GRID* grid() const;
		double energy() const;
	protected:
		void __integrate(double&, double&, T&...);
		GRID *gr;
		double e;
		Func<double, T...> V;
	private:
		int Z, n, l;
		ArrayXd uu;
};

using VHZorbHOM = Verlet_HZorb<HomogeneousGrid, int, int, double>;

class Verlet_HZorb_HOM : public VHZorbHOM{
	public:
		Verlet_HZorb_HOM(int, int, int, int, double, double);
	private:
		static double __V(double&, double&, int&, int&, double&);
};

using VHZorbTHI = Verlet_HZorb<ThijssenGrid, int, int, ThijssenGrid, double>;

class Verlet_HZorb_THI : public VHZorbTHI{
	public:
		Verlet_HZorb_THI(int, int, int, int, double, double, double);
	private:
		static double __V(double&, double&, int&, int&, ThijssenGrid&, double&);
};

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/

Grid::Grid(unsigned long nmax, double rmax):nmax(nmax), rmax(rmax), h(rmax/nmax){}

Grid::Grid(const Grid &o): nmax(o.nmax), rmax(o.rmax), h(o.h), r(o.r){}

HomogeneousGrid::HomogeneousGrid(unsigned long nmax, double rmax): Grid(nmax, rmax){
	Grid::r = ArrayXd::LinSpaced(nmax, 0., nmax)*Grid::h + EPS;
}

HomogeneousGrid::HomogeneousGrid(const HomogeneousGrid &o): Grid(o){}

ArrayXd HomogeneousGrid::jacobian() const{ return ArrayXd::Ones(Grid::nmax); }

ArrayXd HomogeneousGrid::wfFactor() const{ return ArrayXd::Ones(Grid::nmax); }

ThijssenGrid::ThijssenGrid(unsigned long nmax, double rmax, 
		double delta): Grid(nmax, rmax), delta(delta){
	Grid::h = 1.;
	rp = rmax/(exp(delta*nmax) - 1.);
	Grid::r = rp*(exp(delta*ArrayXd::LinSpaced(nmax, 0., nmax-1)) - 1. + EPS);
}

ThijssenGrid::ThijssenGrid(const ThijssenGrid &o): Grid(o), delta(o.delta), rp(o.rp){}

ArrayXd ThijssenGrid::jacobian() const{ 
	return rp*delta*exp(delta*ArrayXd::LinSpaced(Grid::nmax, 0., Grid::nmax-1));
}

ArrayXd ThijssenGrid::wfFactor() const{
	return exp(.5*delta*ArrayXd::LinSpaced(Grid::nmax, 0., Grid::nmax-1));
}

template <typename GRID, class ...T>
Verlet_HZorb<GRID, T...>::Verlet_HZorb(int Z, int n, int l): Z(Z), n(n), l(l){}

template <typename GRID, class ...T>
Verlet_HZorb<GRID, T...>::Verlet_HZorb(const Verlet_HZorb &o):
	Z(o.Z), n(o.n), l(o.l), V(*o.V), e(o.e){}

template <typename GRID, class ...T>
Verlet_HZorb<GRID, T...>::~Verlet_HZorb(){delete gr;}

template <typename GRID, class ...T> void Verlet_HZorb<GRID, T...>::__integrate(double &step, 
	double &accuracy, T& ...par){
	double emin=-2.*pow(Z, 2), emax=EPS;
	ArrayXd uuu, j = ArrayXd::LinSpaced(gr->nmax, 0, gr->nmax-1);
	while (abs(emax-emin)>accuracy){
		e = .5*(emax+emin);
		uu = Verlet<double, T...>(0., 1e-10, gr->r, gr->h, Mode::BACKWARD, V, par...);
		uu *= gr->wfFactor();
		int nodes = 0;
		for (int i=0; i<gr->nmax-1; ++i)
			if (uu[i]*uu[i+1] < 0.) nodes += 1;
		if (nodes > n-l-1) emax = e;
		else emin = e;
	}
	uu /= sqrt(simpson(gr->jacobian()*pow(uu, 2), gr->h));
}


template <typename GRID, class ...T> ArrayXd Verlet_HZorb<GRID, T...>::r() const{ return gr->r; }

template <typename GRID, class ...T> ArrayXd Verlet_HZorb<GRID, T...>::u() const{ return uu; }

template <typename GRID, class ...T> GRID* Verlet_HZorb<GRID, T...>::grid() const{ 
	return new GRID(*gr);
}

template <typename GRID, class ...T> double Verlet_HZorb<GRID, T...>::energy() const{ return e; }

Verlet_HZorb_HOM::Verlet_HZorb_HOM(int Z, int n, int l, int nmax, 
		double rmax, double accuracy):VHZorbHOM(Z, n, l){
	VHZorbHOM::gr = new HomogeneousGrid(nmax, rmax);
	VHZorbHOM::V = __V;
	VHZorbHOM::__integrate(VHZorbHOM::gr->h, accuracy, Z, l, VHZorbHOM::e);
	}

double Verlet_HZorb_HOM::__V(double &r, double &y, int &ZZ, int &ll, double &ee){
	return 2.*((.5*ll*(ll+1) - r*ZZ)/pow(r, 2) - ee)*y;
}



Verlet_HZorb_THI::Verlet_HZorb_THI(int Z, int n, int l, int nmax, 
		double rmax, double delta, double accuracy):VHZorbTHI(Z, n, l){
		VHZorbTHI::gr = new ThijssenGrid(nmax, rmax, delta);
		VHZorbTHI::V = __V;
		VHZorbTHI::__integrate(VHZorbTHI::gr->h, accuracy, 
				Z, l, *VHZorbTHI::gr, VHZorbTHI::e);
	}

double Verlet_HZorb_THI::__V(double &r, double &y, int &ZZ, int &ll, 
		ThijssenGrid &ggr, double &ee){
	return pow(ggr.delta, 2)*(.25 + 2.*pow(r+ggr.rp, 2)*((.5*ll*(ll+1) 
					- r*ZZ)/pow(r, 2) - ee))*y;
}

#endif
