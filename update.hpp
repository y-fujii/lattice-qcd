#pragma once

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

	/*
	return MatrixSU2<double>(
		cos( th0 ),
		sin( th0 ) * cos( th1 ),
		sin( th0 ) * sin( th1 ) * cos( th2 ),
		sin( th0 ) * sin( th1 ) * sin( th2 )
	);
	*/

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

struct ApxMutator {
	ApxMutator( double e ): _eps( e ) {}

	template<class RandGen>
	MatrixSU2<double> operator()( MatrixSU2<double> const& u, RandGen& randGen ) {
		int axis = uniform_int<int>( 0, 2 )( randGen );
		double t = uniform_real_fast<double>( -_eps, +_eps )( randGen );
		MatrixSU2<double> m( sqrt( 1.0 - t * t ), 0.0, 0.0, 0.0 );
		*(&m.a1 + axis) = t;

		return m * u;
	}

	private:
		double const _eps;
};

template<class Lattice, class Mutator, class RandGen>
inline bool updateLocal( Lattice& lat, Site<Lattice::ndim> const& x, int mu, double beta, Mutator mutate, RandGen& randGen ) {
	typename Lattice::Elem v = Lattice::Elem::zero();
	for( int nu = 0; nu < Lattice::ndim; ++nu ) {
		if( mu == nu ) continue;
		v = v + lat( x + mu, nu ) * inv( lat( x, nu ) * lat( x + nu, mu ) );
		v = v + inv( lat( x - nu, mu ) * lat( x + mu - nu, nu ) ) * lat( x - nu, nu );
	}

	typename Lattice::Elem u0 = lat( x, mu );
	typename Lattice::Elem u1 = mutate( u0, randGen );
	double dE = tr( (u1 - u0) * v );
	double p = uniform_real_fast<double>( 0.0, 1.0 )( randGen );
	if( p <= exp( beta * dE ) ) {
		lat( x, mu ) = u1;
		return true;
	}
	else {
		return false;
	}
}

template<class Lattice, class Mutator, class RandGen>
inline bool update0( Lattice& lat, double beta, Mutator mutate, RandGen& randGen ) {
	Site<Lattice::ndim> x;
	for( int i = 0; i < Lattice::ndim; ++i ) {
		x[i] = uniform_int<int>( 0, lat.size() - 1 )( randGen );
	}
	int mu = uniform_int<int>( 0, Lattice::ndim - 1 )( randGen );

	return updateLocal( lat, x, mu, beta, mutate, randGen );
}

template<class Lattice, class Mutator, class RandGen>
inline double update1( Lattice& lat, double beta, Mutator mutate, RandGen& randGen ) {
	int nAcc = 0;
	Site<Lattice::ndim> x;
	x.assign( 0 );
	while( lat.next( x ) ) {
		for( int mu = 0; mu < Lattice::ndim; ++mu ) {
			if( updateLocal( lat, x, mu, beta, mutate, randGen ) ) {
				++nAcc;
			}
		}
	}
	return double( nAcc ) / lat.nLinks();
}
