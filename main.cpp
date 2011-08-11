#include <algorithm>
#include <iostream>
#include <tr1/functional>
#include <tr1/random>
#include <fenv.h>

#include <Eigen/Core>
#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION( Eigen::Matrix2cd )
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION( Eigen::Matrix3cd )

#include "matrix-su2.hpp"
#include "matrix-sun.hpp"
#include "mutator.hpp"
#include "link-lattice.hpp"
#include "update.hpp"
#include "wilson.hpp"

using namespace std;


int main() {
	feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW );

	int const L = 8;
	double const beta = 5.5;
	LinkLattice<Matrix<complex<double>, 3, 3>, 4> lat( L * 2 );
	tr1::mt19937 rng;

	//generate( lat.begin(), lat.end(), bind( randSU2<mt19937>, 1.0, rng ) );
	//generate( lat.begin(), lat.end(), MatrixSU2<double>::one );
	fill( lat.begin(), lat.end(), one( lat( 0 ) ) );

	for( int i = 0; i < 32; ++i ) {
		double p = 0.0;
		for( int j = 0; j < 16; ++j ) {
			p += update( lat, GaugeAction( beta ), SUnMutator( 0.3 ), rng );
		}
		cerr
			<< get<0>( avgWilsonLoop( lat, 1, 1 ) ) << " "
			<< p / 16 << endl;
	}

	int    wcnt[L][L];
	double wval[L][L];
	double wvar[L][L];
	fill( &wcnt[0][0], &wcnt[L][L], 0 );
	fill( &wval[0][0], &wval[L][L], 0.0 );
	fill( &wvar[0][0], &wvar[L][L], 0.0 );
	for( int i = 0; i < 1024; ++i ) {
		for( int j = 0; j < 64; ++j ) {
			update( lat, GaugeAction( beta ), SUnMutator( 0.3 ), rng );
		}
		for( int j = 1; j <= L; ++j ) {
			for( int k = 1; k < j; ++k ) {
				double s01, s02;
				double s11, s12;
				tie( s01, s02 ) = avgWilsonLoop( lat, j, k );
				tie( s11, s12 ) = avgWilsonLoop( lat, k, j );
				wval[j-1][k-1] += s01 + s11;
				wvar[j-1][k-1] += s02 + s12;
				wcnt[j-1][k-1] += 2;
			}
			double s1, s2;
			tie( s1, s2 ) = avgWilsonLoop( lat, j, j );
			wval[j-1][j-1] += s1;
			wvar[j-1][j-1] += s2;
			wcnt[j-1][j-1] += 1;
		}
		cerr << get<0>( avgWilsonLoop( lat, 1, 1 ) ) << endl;
	}

	for( int j = 1; j <= L; ++j ) {
		for( int k = 1; k < j; ++k ) {
			double w  = wval[j-1][k-1] / wcnt[j-1][k-1];
			double dw = sqrt( (wvar[j-1][k-1] / wcnt[j-1][k-1] - w * w) / (wcnt[j-1][k-1] * lat.nSites() * 6 - 1) );
			cout << j * k << " " << w << " " << dw << "\n";
		}
		double w  = wval[j-1][j-1] / wcnt[j-1][j-1];
		double dw = sqrt( (wvar[j-1][j-1] / wcnt[j-1][j-1] - w * w) / (wcnt[j-1][j-1] * lat.nSites() * 6 - 1) );
		cout << j * j << " " << w << " " << dw << "\n";
	}
	cout << endl;

	return 0;
}
