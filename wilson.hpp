#pragma once

#include "matrix-su2.hpp"
#include "link-lattice.hpp"


template<class Matrix>
inline double plaquette( LinkLattice<Matrix>& lat, Site<> const& x, int mu, int nu ) {
	return tr(
		   ( lat( x, mu ) * lat( x + mu, nu ) ) *
		inv( lat( x, nu ) * lat( x + nu, mu ) )
	);
}

template<class Matrix>
inline double wilsonLoop( LinkLattice<Matrix>& lat, Site<> const& x, int mu, int nu, int w, int h ) {
	Site<> y = x;
	Matrix u = Matrix::one();
	for( int i = 0; i < w; ++i ) {
		u = u * lat( y, mu );
		y = y + mu;
	}
	for( int i = 0; i < h; ++i ) {
		u = u * lat( y, nu );
		y = y + nu;
	}

	Site<> z = x;
	Matrix v = Matrix::one();
	for( int i = 0; i < h; ++i ) {
		v = v * lat( z, nu );
		z = z + nu;
	}
	for( int i = 0; i < w; ++i ) {
		v = v * lat( z, mu );
		z = z + mu;
	}
			
	return tr( u * inv( v ) );
}

template<class Matrix>
double avgWilsonLoop( LinkLattice<Matrix>& lat, int w, int h ) {
	double s = 0.0;
	Site<> x = Site<>::zero();
	while( lat.next( x ) ) {
		for( int mu = 0; mu < 4; ++mu ) {
			for( int nu = 0; nu < mu; ++nu ) { 
				s += wilsonLoop( lat, x, mu, nu, w, h );
			}
		}
	}

	return s / (lat.nSites() * 6);
}
