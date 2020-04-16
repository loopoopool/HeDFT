#include<iostream>
#include "helium.h"
using namespace std;

int main(){
	int n;
	double rmax, acc, delta;
	cout << "\n**************************************************\n\n";
	cout << "Insert the following parameters:\n";
	cout << "--> Number of grid points nmax > ", cin >> n;
	cout << "--> Maximum radial distance rmax> ", cin >> rmax;
	cout << "--> Accuracy > ", cin >> acc;
	cout << "--> Parameter delta > ", cin >> delta;
	HeDFT_THI he(n, rmax, acc, delta);
	cout << "\nKohn-Sham eigenvalue = " << he.KSenergy() << endl;
	cout << "Energy = " << he.energy() << endl;
	cout << "\n**************************************************\n\n";
	return 0;
}
