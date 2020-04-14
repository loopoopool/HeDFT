#ifndef ATOM_H
#define ATOM_H

#include "hydrogenic.h"
#include <cmath>

/*************************************************
 * CLASSES DECLARATIONS
 *************************************************/

template <typename GRID> class Hartree{
	public:
		Hartree() = default;
		Hartree(double, double, GRID);
		Hartree(const Hartree&);
		virtual ~Hartree() = default;
		ArrayXd operator()(ArrayXd&);
		ArrayXd operator()() const;
	protected:
		Func<int, ArrayXd, GRID> f;
	private:
		int nmax;
		double rmax, Urmax, accuracy;
		GRID gr;
		ArrayXd VH;
};

class Hartree_HOM : public Hartree<HomogeneousGrid>{
	public:
		Hartree_HOM(double, double, HomogeneousGrid&);
		Hartree_HOM(const Hartree_HOM&) = default;
		~Hartree_HOM() = default;
	private:
		static double __HartreeFunc(int&, double&, 
				ArrayXd&, HomogeneousGrid&);
};

class Hartree_THI : public Hartree<ThijssenGrid>{
	public:
		Hartree_THI(double, double, ThijssenGrid&);
		Hartree_THI(const Hartree_THI&) = default;
		~Hartree_THI() = default;
	private:
		static double __HartreeFunc(int&, double&, 
				ArrayXd&, ThijssenGrid&);
};

template<typename GRID, typename ...T> class HeDFT{
	public:
		HeDFT():E(EPS){};
		HeDFT(const HeDFT&);
		virtual ~HeDFT();
		double energy() const;
	protected:
		ArrayXd __vxc_CeperleyAdler() const;
		void __KS_solve(double&);
		const unsigned int Z=2;
		double E;
		const double CA_A_u{.0311}, CA_B_u{-.048}, CA_C_u{.002}, CA_D_u{.0116},
		      CA_gamma_u{-.1423}, CA_beta1_u{1.0529}, CA_beta2_u{.3334};
		ArrayXd ns;
		Verlet_HZorb<GRID, T...> *hzo;
		Hartree<GRID> *V_H;
		Func<int, ArrayXd, GRID, double> f;
};

using HeDFTHOM = HeDFT<HomogeneousGrid, int, int, double>;

struct HeDFT_HOM : public HeDFTHOM{
	HeDFT_HOM(int, double, double);
	static double __HeKSfunc(int&, double&, ArrayXd&, HomogeneousGrid&, double&);
};

using HeDFTTHI = HeDFT<ThijssenGrid, int, int, ThijssenGrid, double>;

struct HeDFT_THI : public HeDFTTHI{
	HeDFT_THI(int, double, double, double);
	static double __HeKSfunc(int&, double&, ArrayXd&, ThijssenGrid&, double&);
};

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/

template <typename GRID>
Hartree<GRID>::Hartree(double Urmax, double aa, GRID gr): 
	Urmax(Urmax), accuracy(aa), gr(gr){}

template <typename GRID>
Hartree<GRID>::Hartree(const Hartree &o): accuracy(o.accuracy), f(o.f), gr(o.gr){}

template <typename GRID>
ArrayXd Hartree<GRID>::operator()(ArrayXd &n){
	double a, amin = 0.,
	       amax = Urmax/gr.rmax + .5*gr.rmax*(pow(gr.r, 2)*n).maxCoeff();
	ArrayXd U(gr.nmax);
	ArrayXi j = ArrayXi::LinSpaced(gr.nmax, 0, gr.nmax-1);
	while (abs(amax-amin)/(amax+amin)*.5 > accuracy){
		a = .5*(amax+amin);
		U = Verlet(EPS, 1e-10, j, gr.h, Mode::FORWARD, f, n, gr);
		U *= gr.wfFactor();
		U += a*gr.r;
		if (U[gr.nmax-1] > Urmax) amax = a;
		else amin = a;
	} 
	VH = U/gr.r;
	return VH;
}

template <typename GRID>
ArrayXd Hartree<GRID>::operator()() const{ return VH; }

Hartree_HOM::Hartree_HOM(double Urmax, double accuracy, HomogeneousGrid &gr):
	Hartree(Urmax, accuracy, gr){ f = __HartreeFunc; }

double Hartree_HOM::__HartreeFunc(int &i, double &U, 
		ArrayXd &n, HomogeneousGrid &gr){ return -gr.r[i]*n[i]; }

Hartree_THI::Hartree_THI(double Urmax, double accuracy, ThijssenGrid &gr):
	Hartree(Urmax, accuracy, gr){ f = __HartreeFunc; }

double Hartree_THI::__HartreeFunc(int &i, double &U, 
		ArrayXd &n, ThijssenGrid &gr){
	return pow(gr.delta, 2)*(.25*U - sqrt(gr.rp)*pow(gr.r[i] + 
				gr.rp, 1.5)*n[i]*gr.r[i]);
}

template <typename GRID, typename ...T>
HeDFT<GRID, T...>::HeDFT(const HeDFT &o): E(o.E), ns(o.ns){
	hzo = new Verlet_HZorb<GRID, T...>(*o.hzo);
	V_H = new Hartree<GRID>(*o.V_H);
}

template <typename GRID, typename ...T>
HeDFT<GRID, T...>::~HeDFT(){ delete hzo; delete V_H; }

template <typename GRID, typename ...T>
ArrayXd HeDFT<GRID, T...>::__vxc_CeperleyAdler() const{
	int nmax = hzo->grid()->nmax;
	double ec;
	const double a13{1./3.}, a23{2./3.}, a76{7./6.}, a43{4./3.};
	ArrayXd Vx = -pow(ns/a43/(M_PI*M_PI), 1./3.);
//		rs = pow(a43*M_PI*nn, -1./3.),
//		Vc(nmax);
//	for (int i=0; i<nmax; ++i)
//		if (abs(rs[i]-1.) < EPS)
//			Vc[i] = CA_A_u*(log(rs[i]) - a13) + CA_B_u + a23*CA_C_u*rs[i]*log(rs[i]) 
//				+ a13*(2.*CA_D_u - CA_C_u)*rs[i];
//		else{
//			ec = CA_gamma_u/(1. + CA_beta1_u*pow(rs[i], 2) + CA_beta2_u*rs[i]);
//			Vc[i] = pow(ec, 2)/CA_gamma_u*(1. + a76*CA_beta1_u*pow(rs[i], 2) 
//					+ a43*CA_beta2_u*rs[i]);
//		}
	return Vx;
}

template <typename GRID, typename ...T>
double HeDFT<GRID, T...>::energy() const{ return E; }

template <typename GRID>
double average(ArrayXd &n, ArrayXd &y, GRID &gr){
	return simpson(gr.jacobian()*pow(gr.r, 2)*n*y, gr.h);
}

template <typename GRID>
double average(ArrayXd &n, ArrayXd &&y, GRID &gr){
	return simpson(gr.jacobian()*pow(gr.r, 2)*n*y, gr.h);
}

template <typename GRID, typename ...T>
void HeDFT<GRID, T...>::__KS_solve(double &accuracy){
	double e{hzo->energy()}, emin, emax{EPS}, E0{EPS};
	int nn = hzo->grid()->nmax;
	unsigned int nodes;
	ArrayXd phi(nn), V(nn), Vi(nn), Vx(nn);
	ArrayXi j = ArrayXi::LinSpaced(nn, 0, nn-1);
	Vx = __vxc_CeperleyAdler();
	V = (*V_H)(ns) + Vx;
	E = 2*e - .5*average(ns, V + .5*Vx, *hzo->grid());
	clog << "E0 = " << E << endl;
	emin = -4.;
	while (abs(e-E0) > accuracy){
		E0 = e;
		e = .5*(emax+emin);
		phi = Verlet(EPS, 1e-10, j, hzo->grid()->h, Mode::BACKWARD, f, V, *hzo->grid(), e);
		phi *= hzo->grid()->wfFactor();
		nodes = 0;
		for (int i=0; i<nn-1; ++i)
			if (phi[i]*phi[i+1] < 0.) nodes += 1;
		if (nodes > 0) emax = e;
		else emin = e;
		phi /= sqrt(simpson(hzo->grid()->jacobian()*pow(phi, 2), hzo->grid()->h));
		ns = 2.*pow(phi/hzo->r(), 2);
		Vx = __vxc_CeperleyAdler();
		V = (*V_H)(ns) + Vx;
		E = 2.*e - .5*average(ns, (*V_H)() + .5*Vx, *hzo->grid());
	}
	clog << "e = " << e << endl;
}

HeDFT_HOM::HeDFT_HOM(int nmax, double rmax, double accuracy){
	hzo = new Verlet_HZorb_HOM(Z, 1, 0, nmax, rmax, accuracy);
	ns = 2.*pow(hzo->u()/hzo->r(), 2);
	V_H = new Hartree_HOM(Z, accuracy, *hzo->grid());
	f = __HeKSfunc;
	__KS_solve(accuracy);
}

double HeDFT_HOM::__HeKSfunc(int &i, double &y, ArrayXd &V, HomogeneousGrid &gr, double &e){
	return 2.*(-2./gr.r[i] + V[i] - e)*y; 
}

HeDFT_THI::HeDFT_THI(int nmax, double rmax, double accuracy, double delta){
	hzo = new Verlet_HZorb_THI(Z, 1, 0, nmax, rmax, accuracy, delta);
	ns = pow(hzo->u()/hzo->r(), 2);
	V_H = new Hartree_THI(Z, accuracy, *hzo->grid());
	f = __HeKSfunc;
	__KS_solve(accuracy);
}

double HeDFT_THI::__HeKSfunc(int &i, double &y, ArrayXd &V, ThijssenGrid &gr, double &e){
	return pow(gr.delta, 2)*(.25 + 2.*pow(gr.r[i] + gr.rp, 2)*(-2./gr.r[i] + V[i] - e))*y;
}

#endif
