#include <algorithm>
#include <tr1/functional>
#include <fenv.h>
#include "matrix-su2.hpp"
#include "link-lattice.hpp"
#include "update.hpp"
#include "wilson.hpp"

using namespace std;
using namespace tr1;


int main() {
	feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW );

	double const beta = 3.0 / 2.0;
	LinkLattice<MatrixSU2<double>, 4> lat( 8 );
	mt19937 randGen;

	generate( lat.begin(), lat.end(), bind( randSU2<mt19937>, 1.0, randGen ) );
	//generate( lat.begin(), lat.end(), MatrixSU2<double>::one );

	for( int i = 0; i < 1024; ++i ) {
		update1( lat, beta, SphMutator( 0.25 ), randGen );
	}

	double w[4][4];
	fill( &w[0][0], &w[0][0] + 4 * 4, 0.0 );
	for( int i = 0; i< 1024 / 16; ++i ) {
		/*
		for( int j = 0; j < lat.nLinks() * 16; ++j ) {
			update0( lat, beta, SphMutator( 0.25 ), randGen );
		}
		*/
		for( int j = 0; j < 16; ++j ) {
			update1( lat, beta, SphMutator( 0.25 ), randGen );
		}
		for( int j = 1; j < 4; ++j )
		for( int k = 1; k <= j; ++k ) {
			w[j-1][k-1] += avgWilsonLoop( lat, j, k );
			w[j-1][k-1] += avgWilsonLoop( lat, k, j );
			cout << j * k << " " << log( w[j-1][k-1] / ((i + 1) * 4) ) << "\n";
		}
		cout << endl;
	}

	return 0;
}
