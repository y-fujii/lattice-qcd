#include <tr1/functional>
#include <fenv.h>
#include "format.hpp"
#include "matrix-su2.hpp"
#include "link-lattice.hpp"
#include "update.hpp"
#include "wilson.hpp"

using namespace tr1;


template<class Lattice, class RandGen>
void init( Lattice& lat, RandGen& randGen ) {
	for( int i = 0; i < lat.nLinks(); ++i ) {
		//lat( i ) = randSU2( 1.0, randGen );
		lat( i ) = Lattice::Elem::one();
	}
}

int main() {
	feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW );

	//double const beta = 3.0 / 2.0;
	double const beta = 2.3 / 2.0;
	LinkLattice<MatrixSU2<double>, 4> lat( 8 );

	mt19937 randGen;
	init( lat, randGen );
	while( true ) {
		/*
		int rate = 0;
		for( int i = 0; i < lat.nLinks(); ++i ) {
			rate += update0( lat, beta, ApxMutator( 0.75 ), randGen ) ? 1 : 0;
			//rate += update0( lat, beta, SphMutator( 0.25 ), randGen ) ? 1 : 0;
		}
		cout << double( rate ) / lat.nLinks() << endl;
		/**/
		cout << update1( lat, beta, ApxMutator( 1.0 ), randGen ) << ": ";
		//update1( lat, beta, SphMutator( 0.25 ), randGen );
		cout
			<< avgWilsonLoop( lat, 1, 1 ) / 2.0 << " "
			<< avgWilsonLoop( lat, 2, 2 ) / 2.0 << " "
			<< avgWilsonLoop( lat, 3, 3 ) / 2.0 << " "
			<< avgWilsonLoop( lat, 4, 4 ) / 2.0 << endl
		;
	}

	return 0;
}
