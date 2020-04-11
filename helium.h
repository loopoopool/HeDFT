#ifndef ATOM_H
#define ATOM_H

#include "hydrogenic.h"
#include <cmath>

/*************************************************
 * CLASSES DECLARATIONS
 *************************************************/

template <typename GRID, typename ...T> class Hartree{
	public:
		Hartree(double&, double&, GRID&);
		Hartree(const Hartree&);
		virtual ~Hartree(){ delete f; }
		ArrayXd operator()(ArrayXd&, T...) const;
	protected:
		virtual double __poisson(double&, double&, T...)=0;
	private:
		Func<T...> *f;
		unsigned long nmax;
		double rmax, Urmax, accuracy;
		GRID *gr;
};

class Hartree_HOM : Hartree<HomogeneousGrid>{
	public:	Hartree_HOM(unsigned long, double, double);
	private: static double __poisson(double&, double&, double&);
};


class Hartree_THI : Hartree<HomogeneousGrid>{
	public: Hartree_THI(unsigned long, double, double);
	private: static double __poisson(double&, double&, double&);
};

template<typename GRID, typename ...T> class HeDFT{
	public:
		HeDFT() = default;
		HeDFT(const HeDFT&);
		virtual ~HeDFT();
		double energy() const;
	protected:
		ArrayXd __vxc_CeperleyAdler() const;
		void __KS_solve(ArrayXd&, double&, double&);
		const unsigned int Z=2;
		double E, emax, emin;
		const double CA_A_u{.0311}, CA_B_u{-.048}, CA_C_u{.002}, CA_D_u{.0116},
		      CA_gamma_u{-.1423}, CA_beta1_u{1.0529}, CA_beta2_u{.3334};
		ArrayXd r, n;
		Verlet_HZorb<GRID, T...> *hzo;
		Hartree<GRID, T...> *V_H;
};

using HDFTHOM = HeDFT<HomogeneousGrid, int&, int&, double&>;

struct HeDFT_HOM : HeDFT<HomogeneousGrid, int&, int&, double&>{
	HeDFT_HOM(unsigned long, double, double);
};

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/


template <typename GRID, typename ...T>
Hartree<GRID, T...>::Hartree(double &Urmax, double &aa, GRID &gr): 
	Urmax(Urmax), accuracy(aa), gr(gr){
		f = new Func<T...>(__poisson);
}

template <typename GRID, typename ...T>
Hartree<GRID, T...>::Hartree(const Hartree &o): accuracy(o.accuracy){
	f = new Func<T...>(o.f); 
}

Hartree_HOM::Hartree_HOM(unsigned long &nmax, double &rmax, double &accuracy)

template <typename GRID, typename ...T>
ArrayXd Hartree<GRID, T...>::operator()(ArrayXd &n, T ...par) const{
	double a, amin = 0.,
	       amax = Urmax/gr->rmax + .5*gr->rmax*(pow(gr->r, 2)*n).maxCoeff();
	ArrayXd U(gr->nmax); 
	while (abs(amax-amin)/(amax+amin)*.5 > accuracy){
		a = .5*(amax+amin);
		U = Verlet<T...>(EPS, 1e-10, gr->r, gr->h, Mode::FORWARD, f, par...);
		U *= gr->wfFactor();
		if (U[gr->nmax-1] > Urmax) amax = a;
		else amin = a;
	}
	return U/gr->r;
}

template <typename GRID, typename ...T>
HeDFT<GRID, T...>::HeDFT(const HeDFT &o): E(o.E), r(o.r), n(o.n){
	hzo = new Verlet_HZorb<GRID, T...>(*o.hzo);
	V_H = new Hartree<GRID, T...>(*o.V_H);
}

template <typename GRID, typename ...T>
HeDFT<GRID, T...>::~HeDFT(){ delete hzo; delete V_H; }

template <typename GRID, typename ...T>
ArrayXd HeDFT<GRID, T...>::__vxc_CeperleyAdler() const{
	unsigned long nmax = hzo->gr->nmax;
	double ec;
	const double a13{1./3.}, a23{2./3.}, a76{7./6.}, a43{4./3.};
	ArrayXd Vx = -pow(3./M_PI*n, 1./3.),
		rs = pow(3./4./M_PI/n, 1./3.),
		Vc(n);
	for (unsigned long i=0; i<nmax; ++i)
		if (abs(rs[i]-1.) < EPS)
			Vc[i] = CA_A_u*(log(rs[i]) - a13) + CA_B_u + a23*CA_C_u*rs[i]*log(rs[i]) 
				+ a13*(2.*CA_D_u - CA_C_u)*rs[i];
		else{
			ec = CA_gamma_u/(1. + CA_beta1_u*pow(rs[i], 2) + CA_beta2_u*rs[i]);
			Vc[i] = pow(ec, 2)/CA_gamma_u*(1. + a76*CA_beta1_u*pow(rs[i], 2) 
					+ a43*CA_beta2_u*rs[i]);
		}
	return Vx + Vc;
}

template <typename GRID, typename ...T>
double HeDFT<GRID, T...>::energy() const{ return E; }

//template <typename GRID, typename ...T>
//void HeDFT<GRID, T...>::__KS_solve(ArrayXd &V_KS0, double &Eprev, double &accuracy){
//	double e;
//	while (abs(E-Eprev) > accuracy){
//		e = .5*(emax+emin)
//		phi = Verlet<T...>(EPS, 1e-10, gr->r, gr->h, Mode::FORWARD, f, par...);
//		phi *= gr->wfFactor();
//		if (U[gr->nmax-1] > Urmax) amax = a;
//		else amin = a;

HeDFT_HOM::HeDFT_HOM(unsigned long nmax, double rmax, double accuracy){
	hzo = new Verlet_HZorb_HOM(HDFTHOM::Z, 1, 0, nmax, rmax, accuracy);
	r = hzo->r();
	n = pow(hzo->u()/r, 2);
	V_H = new Hartree_HOM(nmax, rmax, accuracy);



#endif
