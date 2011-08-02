#pragma once

#include <iterator>
#include <iostream>
#include <cmath>
#include <tr1/random>
#include "matrix-su2.hpp"
#include "link-lattice.hpp"

using namespace std;
using namespace tr1;


template<class T>
struct uniform_real_fast {
	uniform_real_fast( T const& min, T const& max ):
		_min( min ), _max( max ) {}

	template<class RandGen>
	T operator()( RandGen& gen ) {
		assert( gen.min() == 0 );
		return gen() * ((_max - _min) / (1.0 + gen.max())) + _min;
	}

	private:
		T const _min, _max;
};

template<class RandGen>
inline MatrixSU2<double> randSU2( double eps, RandGen& randGen ) {
	double th0 = uniform_real_fast<double>( 0.0, M_PI * eps )( randGen );
	double th1 = uniform_real_fast<double>( 0.0, M_PI )( randGen );
	double th2 = uniform_real_fast<double>( 0.0, M_PI * 2.0 )( randGen );

	double c0 = cos( th0 );
	double s0 = sqrt( 1 - c0 * c0 );
	double c1 = cos( th1 );
	double s0s1 = s0 * sqrt( 1 - c1 * c1 );

	return MatrixSU2<double>(
		c0,
		s0 * c1,
		s0s1 * cos( th2 ),
		s0s1 * sin( th2 )
	);
}

struct SphMutator {
	SphMutator( double e ): _eps( e ) {}

	template<class RandGen>
	MatrixSU2<double> operator()( MatrixSU2<double> const& u, RandGen& randGen ) {
		return randSU2( _eps, randGen ) * u;
	}

	private:
		double const _eps;
};

struct GaugeAction {
	GaugeAction( double g ): _g( g ) {}

	template<class Lattice>
	double operator()( Lattice& lat, Site<Lattice::ndim> const& x, int mu, typename Lattice::Elem const& u ) {
		typename Lattice::Elem v = Lattice::Elem::zero();
		for( int nu = 0; nu < Lattice::ndim; ++nu ) {
			if( mu == nu ) continue;
			Site<Lattice::ndim> y = x - nu;
			v = v + lat( x + mu, nu ) * inv( lat( x, nu ) * lat( x + nu, mu ) );
			v = v + inv( lat( y, mu ) * lat( y + mu, nu ) ) * lat( y, nu );
		}

		return (-_g) * tr( (u - lat( x, mu )) * v );
	}

	private:
		double const _g;
};

template<class Lattice, class ActionDiff, class Mutator, class RandGen>
double update( Lattice& lat, ActionDiff actionDiff, Mutator mutate, RandGen& randGen ) {
	int nAcc = 0;
	Site<Lattice::ndim> x;
	x.assign( 0 );
	while( lat.next( x ) ) {
		for( int mu = 0; mu < Lattice::ndim; ++mu ) {
			typename Lattice::Elem u = mutate( lat( x, mu ), randGen );
			double dS = actionDiff( lat, x, mu, u );
			double p = uniform_real_fast<double>( 0.0, 1.0 )( randGen );
			if( p <= exp( -dS ) ) {
				lat( x, mu ) = u;
				++nAcc;
			}
		}
	}
	return double( nAcc ) / lat.nLinks();
}
