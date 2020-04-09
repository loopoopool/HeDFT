#ifndef INTER_H
#define INTER_H

#include <iostream>
#include <cmath>
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

template <class ...T> 
using VerletField = double(*)(double&, double&, T...);

template <class ...T> ArrayXd Verlet(double y0, double y1, ArrayXd &t, 
		double step, Mode m, VerletField<T...> f, T... par){
	int n = t.size();
	ArrayXd tt(n), y(n);
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
	Grid(int, double);
	Grid(int, double, double);
	Grid(const Grid &) = default;
	~Grid() = default;
	ArrayXd r;
	double rp, delta, rmax, h;
	int nmax;
	GridType gt;
};

template <GridType gt, typename ...T> class Verlet_HZorb{
	public:
		Verlet_HZorb(int, int, int);
		Verlet_HZorb(const Verlet_HZorb &);
		~Verlet_HZorb();
		ArrayXd u() const;
		ArrayXd r() const;
		double energy() const;
	protected:
		void __integrate(double&, double&, T...);
		Grid *gr;
		VerletField<T...> f;
		double e;
	private:
		int Z, n, l;
		ArrayXd uu;
};

using VHZorbHOM = Verlet_HZorb<GridType::HOMOGENEOUS, int&, int&, double&>;

class Verlet_HZorb_HOM : public VHZorbHOM{
	public:
		Verlet_HZorb_HOM(int, int, int, int, double, double);
	private:
		static double __Vhom(double&, double&, int&, int&, double&);
};

using VHZorbTHI = Verlet_HZorb<GridType::THIJSSEN, int&, int&, Grid*, double&>;

class Verlet_HZorb_THI : public VHZorbTHI{
	public:
		Verlet_HZorb_THI(int, int, int, int, double, double, double);
	private:
		static double __Vthijssen(double &, double &, int&, int&, Grid*, double&);
};

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/

Grid::Grid(int nmax, double rmax):nmax(nmax), rmax(rmax){
	gt = GridType::HOMOGENEOUS;
	rp = 0.;
	delta = 0.;
	h = rmax/nmax;
	r = ArrayXd::LinSpaced(nmax, 0., nmax-1)*h + EPS;
}

Grid::Grid(int nmax, double rmax, double delta):
nmax(nmax), rmax(rmax), delta(delta){
	gt = GridType::THIJSSEN;
	h = 1.;
	rp = rmax/(exp(delta*nmax) - 1.);
	r = rp*(exp(delta*ArrayXd::LinSpaced(nmax, 0., nmax-1)) - 1. + EPS);
}

template <GridType gt, class ...T>
Verlet_HZorb<gt, T...>::Verlet_HZorb(int Z, int n, int l): Z(Z), n(n), l(l){}

template <GridType gt, class ...T>
Verlet_HZorb<gt, T...>::Verlet_HZorb(const Verlet_HZorb &o):
	Z(o.Z), n(o.n), l(o.l), e(o.e), f(o.f){ gr = new Grid(*o.gr); }

template <GridType gt, class ...T>
Verlet_HZorb<gt, T...>::~Verlet_HZorb(){delete gr;}

template <GridType gt, class ...T> void Verlet_HZorb<gt, T...>::__integrate(double &step, 
	double &accuracy, T... par){
	double emin=-2.*pow(Z, 2), emax=EPS;
	ArrayXd uuu, j = ArrayXd::LinSpaced(gr->nmax, 0, gr->nmax-1);
	while (abs(emax-emin)>accuracy){
		e = .5*(emax+emin);
		uu = Verlet<T...>(0., 1e-10, gr->r, gr->h, Mode::BACKWARD, f, par...);
		if (gt == GridType::THIJSSEN) uu *= exp(.5*j*gr->delta);
		int nodes = 0;
		for (int i=0; i<gr->nmax-1; ++i)
			if (uu[i]*uu[i+1] < 0.) nodes += 1;
		if (nodes > n-l-1) emax = e;
		else emin = e;
	}
	if (gt == GridType::THIJSSEN) 
		uuu = gr->rp*gr->delta*exp(j*gr->delta)*pow(uu, 2);
	else uuu = pow(uu, 2);
	uu /= sqrt(simpson(uuu, gr->h));
}


template <GridType gt, class ...T> ArrayXd Verlet_HZorb<gt, T...>::r() const{ return gr->r; }

template <GridType gt, class ...T> ArrayXd Verlet_HZorb<gt, T...>::u() const{ return uu; }

template <GridType gt, class ...T> double Verlet_HZorb<gt, T...>::energy() const{ return e; }

Verlet_HZorb_HOM::Verlet_HZorb_HOM(int Z, int n, int l, int nmax, 
		double rmax, double accuracy):VHZorbHOM(Z, n, l){
	VHZorbHOM::gr = new Grid(nmax, rmax);
	VHZorbHOM::f = __Vhom;
	VHZorbHOM::__integrate(VHZorbHOM::gr->h, accuracy, Z, l, VHZorbHOM::e);
	}

double Verlet_HZorb_HOM::__Vhom(double &r, double &y, int &ZZ, int &ll, double &ee){
	return 2.*((.5*ll*(ll+1) - r*ZZ)/pow(r, 2) - ee)*y;
}


Verlet_HZorb_THI::Verlet_HZorb_THI(int Z, int n, int l, int nmax, 
		double rmax, double accuracy, double delta):VHZorbTHI(Z, n, l){
		VHZorbTHI::gr = new Grid(nmax, rmax, delta);
		VHZorbTHI::f = __Vthijssen;
		VHZorbTHI::__integrate(VHZorbTHI::gr->h, accuracy, 
				Z, l, VHZorbTHI::gr, VHZorbTHI::e);
	}

double Verlet_HZorb_THI::__Vthijssen(double &r, double &y, int &ZZ, int &ll, 
		Grid *ggr, double &ee){
	return pow(ggr->delta, 2)*(.25 + 2.*pow(r+ggr->rp, 2)*((.5*ll*(ll+1) 
					- r*ZZ)/pow(r, 2) - ee))*y;
}
#endif
