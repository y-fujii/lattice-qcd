#pragma once

#include <complex>
#include <tr1/tuple>
#include "misc.hpp"

using namespace std;
using namespace tr1;


template<class Lattice>
inline double wilsonLoop( Lattice& lat, Site<Lattice::ndim> const& x, int mu, int nu, int w, int h ) {
	Site<Lattice::ndim> y = x;
	typename Lattice::Elem u = one( typename Lattice::Elem() );
	for( int i = 0; i < w; ++i ) {
		u = u * lat( y, mu );
		y = y + mu;
	}
	for( int i = 0; i < h; ++i ) {
		u = u * lat( y, nu );
		y = y + nu;
	}

	Site<Lattice::ndim> z = x;
	typename Lattice::Elem v = one( typename Lattice::Elem() );
	for( int i = 0; i < h; ++i ) {
		v = v * lat( z, nu );
		z = z + nu;
	}
	for( int i = 0; i < w; ++i ) {
		v = v * lat( z, mu );
		z = z + mu;
	}
			
	return real( ntrace( u * inv( v ) ) );
}

template<class Lattice>
tuple<double, double> avgWilsonLoop( Lattice& lat, int w, int h ) {
	double s1 = 0.0;
	double s2 = 0.0;
	Site<Lattice::ndim> x;
	x.assign( 0 );
	do {
		for( int mu = 0; mu < Lattice::ndim; ++mu ) {
			for( int nu = 0; nu < mu; ++nu ) { 
				double v = wilsonLoop( lat, x, mu, nu, w, h );
				s1 += v;
				s2 += v * v;
			}
		}
	} while( lat.next( x ) );

	int const n = lat.nSites() * (Lattice::ndim * (Lattice::ndim - 1) / 2);
	return make_tuple( s1 / n, s2 / n );
}
