#ifndef GRID_H
#define GRID_H
#include "core.h"

/*************************************************
 * CLASSES DECLARATIONS
 *************************************************/

enum GridType{ HOMOGENEOUS, THIJSSEN };

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

/*************************************************
 * METHODS DEFINITIONS
 *************************************************/

Grid::Grid(unsigned long nmax, double rmax):nmax(nmax), rmax(rmax), 
	h(rmax/nmax){}

Grid::Grid(const Grid &o): nmax(o.nmax), rmax(o.rmax), h(o.h), r(o.r){}

HomogeneousGrid::HomogeneousGrid(unsigned long nmax, double rmax): 
	Grid(nmax, rmax){
		Grid::r = ArrayXd::LinSpaced(nmax, EPS, nmax-1)*Grid::h + EPS;
}

HomogeneousGrid::HomogeneousGrid(const HomogeneousGrid &o): Grid(o){}

ArrayXd HomogeneousGrid::jacobian() const{ return ArrayXd::Ones(Grid::nmax); }

ArrayXd HomogeneousGrid::wfFactor() const{ return ArrayXd::Ones(Grid::nmax); }

ThijssenGrid::ThijssenGrid(unsigned long nmax, double rmax, 
		double delta): Grid(nmax, rmax), delta(delta){
	Grid::h = 1.;
	rp = rmax/(exp(delta*nmax) - 1.);
	Grid::r = rp*(exp(delta*ArrayXd::LinSpaced(nmax, 0., 
					nmax-1)) - 1. + EPS);
}

ThijssenGrid::ThijssenGrid(const ThijssenGrid &o): Grid(o), delta(o.delta), 
	rp(o.rp){}

ArrayXd ThijssenGrid::jacobian() const{ 
	return rp*delta*exp(delta*ArrayXd::LinSpaced(Grid::nmax, 0., 
				Grid::nmax-1));
}

ArrayXd ThijssenGrid::wfFactor() const{
	return exp(.5*delta*ArrayXd::LinSpaced(Grid::nmax, 0., Grid::nmax-1));
}

#endif
