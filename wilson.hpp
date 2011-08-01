#pragma once

#include "matrix-su2.hpp"
#include "link-lattice.hpp"


template<class Lattice>
inline double plaquette( Lattice& lat, Site<Lattice::ndim> const& x, int mu, int nu ) {
	return tr(
		   ( lat( x, mu ) * lat( x + mu, nu ) ) *
		inv( lat( x, nu ) * lat( x + nu, mu ) )
	);
}

template<class Lattice>
inline double wilsonLoop( Lattice& lat, Site<Lattice::ndim> const& x, int mu, int nu, int w, int h ) {
	Site<Lattice::ndim> y = x;
	typename Lattice::Elem u = Lattice::Elem::one();
	for( int i = 0; i < w; ++i ) {
		u = u * lat( y, mu );
		y = y + mu;
	}
	for( int i = 0; i < h; ++i ) {
		u = u * lat( y, nu );
		y = y + nu;
	}

	Site<Lattice::ndim> z = x;
	typename Lattice::Elem v = Lattice::Elem::one();
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

template<class Lattice>
double avgWilsonLoop( Lattice& lat, int w, int h ) {
	double s = 0.0;
	Site<Lattice::ndim> x;
	x.assign( 0 );
	while( lat.next( x ) ) {
		for( int mu = 0; mu < Lattice::ndim; ++mu ) {
			for( int nu = 0; nu < mu; ++nu ) { 
				s += wilsonLoop( lat, x, mu, nu, w, h );
			}
		}
	}

	int const n = Lattice::ndim * (Lattice::ndim - 1) / 2;
	return s / (lat.nSites() * n);
}
