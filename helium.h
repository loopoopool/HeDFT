#ifndef ATOM_H
#define ATOM_H
#include "hydrogenic.h"
#include "hartree.h"

namespace numc{
	const double a13{1./3.}, a23{2./3.}, a76{7./6.}, a43{4./3.},
	      CA_A_u{.0311}, CA_B_u{-.048}, CA_C_u{.002}, CA_D_u{.0116},
	      CA_gamma_u{-.1423}, CA_beta1_u{1.0529}, CA_beta2_u{.3334},
	      aVx{3./(4.*M_PI*M_PI)};
}

/*************************************************
 * CLASSES DECLARATIONS
 *************************************************/

template<typename GRID, typename ...T> class HeDFT{
	public:
		HeDFT():E(EPS){};
		HeDFT(const HeDFT&);
		virtual ~HeDFT();
		double energy() const;
		double KSenergy() const;
	protected:
		ArrayXd __ec() const;
		ArrayXd __Vc_CA() const;
		ArrayXd __Vx_CA() const;
		void __KS_solve(double&);
		const unsigned int Z=2;
		double E, eKS;
		ArrayXd n;
		Verlet_HZorb<GRID, T...> *hzo;
		Hartree<GRID> *V_H;
		Func<int, ArrayXd, GRID, double> f;
};

using HeDFTHOM = HeDFT<HomogeneousGrid, int, int, double>;

struct HeDFT_HOM : public HeDFTHOM{
	HeDFT_HOM(int, double, double);
	static double __HeKSfunc(int&, double&, ArrayXd&, HomogeneousGrid&, 
			double&);
};

using HeDFTTHI = HeDFT<ThijssenGrid, int, int, ThijssenGrid, double>;

struct HeDFT_THI : public HeDFTTHI{
	HeDFT_THI(int, double, double, double);
	static double __HeKSfunc(int&, double&, ArrayXd&, ThijssenGrid&, 
			double&);
};

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/


template <typename GRID, typename ...T>
HeDFT<GRID, T...>::HeDFT(const HeDFT &o): E(o.E), n(o.n){
	hzo = new Verlet_HZorb<GRID, T...>(*o.hzo);
	V_H = new Hartree<GRID>(*o.V_H);
}

template <typename GRID, typename ...T>
HeDFT<GRID, T...>::~HeDFT(){ delete hzo; delete V_H; }

template <typename GRID, typename ...T>
ArrayXd HeDFT<GRID, T...>::__ec() const{
	int nmax = hzo->grid()->nmax;
	ArrayXd rs = pow(n/3., -numc::a13), ec(nmax);
	for (int i=0; i<nmax; ++i)
		if (abs(rs[i]-1.) < EPS) 
			ec[i] = (numc::CA_A_u + rs[i]*numc::CA_C_u)*log(rs[i]) + 
				numc::CA_B_u + rs[i]*numc::CA_D_u;
		else ec[i] = numc::CA_gamma_u/(1. + numc::CA_beta1_u*pow(rs[i], 2) + 
					numc::CA_beta2_u*rs[i]);
	return ec;
}

template <typename GRID, typename ...T>
ArrayXd HeDFT<GRID, T...>::__Vx_CA() const{ 
	return -pow(numc::aVx*n, numc::a13); 
}

template <typename GRID, typename ...T>
ArrayXd HeDFT<GRID, T...>::__Vc_CA() const{
	int nmax = hzo->grid()->nmax;
	ArrayXd rs = pow(n/3., -numc::a13), Vc(nmax);
	for (int i=0; i<nmax; ++i)
		if (abs(rs[i]-1.) < EPS)
			Vc[i] = numc::CA_A_u*(log(rs[i]) - numc::a13) + 
				numc::CA_B_u + numc::a23*numc::CA_C_u*rs[i] * 
				log(rs[i]) + numc::a13*(2.*numc::CA_D_u - 
						numc::CA_C_u)*rs[i];
		else Vc[i] = numc::CA_gamma_u*(1. + numc::a76*numc::CA_beta1_u *
				pow(rs[i], 2) + numc::a43*numc::CA_beta2_u *
				rs[i])*pow(1. + numc::CA_beta1_u*pow(rs[i], 2) + 
						numc::CA_beta2_u*rs[i], -2);
	return Vc;
}

template <typename GRID, typename ...T>
double HeDFT<GRID, T...>::energy() const{ return E; }

template <typename GRID, typename ...T>
double HeDFT<GRID, T...>::KSenergy() const{ return eKS; }

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
	double emin{-4.}, emax{EPS}, E0{EPS}, 
	       phi1{hzo->grid()->rmax/hzo->grid()->nmax};
	int nn = hzo->grid()->nmax;
	ArrayXd phi(nn), V(nn), Vx(nn), Vc(nn), wf{hzo->grid()->wfFactor()};
	ArrayXi j = ArrayXi::LinSpaced(nn, 0, nn-1);
	eKS = hzo->energy();
	Vx = __Vx_CA();
	Vc = __Vc_CA();
	V = (*V_H)(n) + Vx + Vc;
	E = 2*eKS - average(n, .5*(*V_H)() + .25*Vx + Vc - __ec(), 
			*hzo->grid());
	while (abs(E-E0) > accuracy){
		E0 = E;
		while(abs(emax-emin) > accuracy){
			eKS = .5*(emax+emin);
			phi = VerletREV(EPS, phi1, j, hzo->grid()->h, f, V, 
					*hzo->grid(), eKS)*wf;
			phi /= sqrt(simpson(hzo->grid()->jacobian()*pow(phi, 2), 
						hzo->grid()->h));
			if (phi[0] > EPS) emin = eKS;
			else emax = eKS;
		}
		emax=EPS;
		emin=-4.;
		n = 2.*pow(phi/hzo->r(), 2);
		Vx = __Vx_CA();
		Vc = __Vc_CA();
		V = (*V_H)(n) + Vx + Vc;
		E = 2.*eKS - average(n, .5*(*V_H)() + .25*Vx + Vc - __ec(), 
				*hzo->grid());
	}
}

HeDFT_HOM::HeDFT_HOM(int nmax, double rmax, double accuracy){
	hzo = new Verlet_HZorb_HOM(Z, 1, 0, nmax, rmax, accuracy);
	n = 2.*pow(hzo->u()/hzo->r(), 2);
	V_H = new Hartree_HOM(Z, accuracy, *hzo->grid());
	f = __HeKSfunc;
	__KS_solve(accuracy);
}

double HeDFT_HOM::__HeKSfunc(int &i, double &y, ArrayXd &V, HomogeneousGrid &gr, 
		double &e){
	return 2.*(-2./gr.r[i] + V[i] - e)*y; 
}

HeDFT_THI::HeDFT_THI(int nmax, double rmax, double accuracy, double delta){
	hzo = new Verlet_HZorb_THI(Z, 1, 0, nmax, rmax, accuracy, delta);
	n = 2.*pow(hzo->u()/hzo->r(), 2);
	V_H = new Hartree_THI(Z, accuracy, *hzo->grid());
	f = __HeKSfunc;
	__KS_solve(accuracy);
}

double HeDFT_THI::__HeKSfunc(int &i, double &y, ArrayXd &V, ThijssenGrid &gr, 
		double &e){
	return pow(gr.delta, 2)*(.25 + 2.*pow(gr.r[i] + gr.rp, 2)*(-2./gr.r[i] 
				+ V[i] - e))*y;
}

#endif
