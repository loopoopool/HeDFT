#ifndef HYDROGENIC_H
#define HYDROGENIC_H
#include "grid.h"

/*************************************************
 * CLASSES DECLARATIONS
 *************************************************/

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
	public: Verlet_HZorb_HOM(int, int, int, int, double, double);
	private: static double __V(double&, double&, int&, int&, double&);
};

using VHZorbTHI = Verlet_HZorb<ThijssenGrid, int, int, ThijssenGrid, double>;

class Verlet_HZorb_THI : public VHZorbTHI{
	public: Verlet_HZorb_THI(int, int, int, int, double, double, double);
	private: static double __V(double&, double&, int&, int&, ThijssenGrid&, 
				 double&);
};

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/

template <typename GRID, class ...T>
Verlet_HZorb<GRID, T...>::Verlet_HZorb(int Z, int n, int l): Z(Z), n(n), l(l){}

template <typename GRID, class ...T>
Verlet_HZorb<GRID, T...>::Verlet_HZorb(const Verlet_HZorb &o):
	Z(o.Z), n(o.n), l(o.l), V(*o.V), e(o.e){}

template <typename GRID, class ...T>
Verlet_HZorb<GRID, T...>::~Verlet_HZorb(){delete gr;}

template <typename GRID, class ...T> 
void Verlet_HZorb<GRID, T...>::__integrate(double &step, 
	double &accuracy, T& ...par){
	double emin{2.*pow(Z, 2)}, emax{EPS}, uu1{pow(gr->rmax/gr->nmax, l+1)};
	ArrayXd j = ArrayXd::LinSpaced(gr->nmax, 0, gr->nmax-1),
		wf = gr->wfFactor();
	while (abs(emax-emin)>accuracy){
		e = .5*(emax+emin);
		uu = VerletREV(EPS, uu1, gr->r, gr->h, V, par...)*wf;
		uu /= sqrt(simpson(gr->jacobian()*pow(uu, 2), gr->h));
		if (uu[0] > EPS) emin = e;
		else emax = e;
	}
}


template <typename GRID, class ...T> 
ArrayXd Verlet_HZorb<GRID, T...>::r() const{ return gr->r; }

template <typename GRID, class ...T> 
ArrayXd Verlet_HZorb<GRID, T...>::u() const{ return uu; }

template <typename GRID, class ...T> 
GRID* Verlet_HZorb<GRID, T...>::grid() const{ return new GRID(*gr); }

template <typename GRID, class ...T> 
double Verlet_HZorb<GRID, T...>::energy() const{ return e; }

Verlet_HZorb_HOM::Verlet_HZorb_HOM(int Z, int n, int l, int nmax, 
		double rmax, double accuracy):VHZorbHOM(Z, n, l){
	VHZorbHOM::gr = new HomogeneousGrid(nmax, rmax);
	VHZorbHOM::V = __V;
	VHZorbHOM::__integrate(VHZorbHOM::gr->h, accuracy, Z, l, VHZorbHOM::e);
}

double Verlet_HZorb_HOM::__V(double &r, double &y, int &ZZ, int &ll, 
		double &ee){
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
