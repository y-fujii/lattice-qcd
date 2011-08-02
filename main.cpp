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
	tr1::mt19937 randGen;

	//generate( lat.begin(), lat.end(), [&]() -> { return randSU2( 1.0, randGen ) } );
	generate( lat.begin(), lat.end(), MatrixSU2<double>::one );

	for( int i = 0; i < 1024; ++i ) {
		update( lat, GaugeAction( beta ), SphMutator( 0.25 ), randGen );
	}

	int    wcnt[16];
	double wval[16];
	double wvar[16];
	fill( wcnt, wcnt + 16, 0 );
	fill( wval, wval + 16, 0.0 );
	fill( wvar, wvar + 16, 0.0 );
	for( int i = 0; i< 512 / 16; ++i ) {
		for( int j = 0; j < 16; ++j ) {
			update( lat, GaugeAction( beta ), SphMutator( 0.25 ), randGen );
		}
		for( int j = 1; j < 4 + 1; ++j ) {
			for( int k = 1; k < 4 + 1; ++k ) {
				double w = avgWilsonLoop( lat, j, k );
				wval[j * k - 1] += w;
				wvar[j * k - 1] += w * w;
				wcnt[j * k - 1] += 1;
			}
		}
		for( int j = 0; j < 16; ++j ) {
			if( wcnt[j] > 0 ) {
				double w  = wval[j] / wcnt[j];
				double dw = sqrt( wvar[j] / wcnt[j] - w * w );
				cout << j + 1 << " " << log( w + dw ) << "\n";
				cout << j + 1 << " " << log( w - dw ) << "\n";
			}
		}
		cout << endl;
	}

	return 0;
}
