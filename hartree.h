#ifndef HARTREE_H
#define HARTREE_H
#include <Eigen/Dense>
#include "grid.h"

using namespace Eigen;

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
	double a, amin, amax, U1{gr.rmax/gr.nmax};
	ArrayXi j = ArrayXi::LinSpaced(gr.nmax, 0, gr.nmax-1);
	ArrayXd U(gr.nmax), wf{gr.wfFactor()};
	amax = Urmax/gr.rmax + gr.rmax*(gr.r*n).maxCoeff();
	amin = Urmax/gr.rmax;
	while (abs(amax-amin) > accuracy){
		a = .5*(amax+amin);
		U = Verlet(EPS, U1, j, gr.h, f, n, gr)*wf;
		U += a*gr.r;
		if (U[gr.nmax-1] - Urmax > accuracy) amax = a;
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

#endif
