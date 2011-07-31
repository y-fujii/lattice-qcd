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
	*/
	return MatrixSU2<double>(
		cos( th0 ),
		sin( th0 ) * cos( th1 ),
		sin( th0 ) * sin( th1 ) * cos( th2 ),
		sin( th0 ) * sin( th1 ) * sin( th2 )
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

template<class Matrix, class Mutator, class RandGen>
inline bool updateLocal( LinkLattice<Matrix>& lat, Site const& x, int mu, double beta, Mutator mutate, RandGen& randGen ) {
	Matrix v = Matrix::zero();
	for( int nu = 0; nu < 4; ++nu ) {
		if( mu == nu ) continue;
		v = v + lat( x + mu, nu ) * inv( lat( x, nu ) * lat( x + nu, mu ) );
		v = v + inv( lat( x - nu, mu ) * lat( x + mu - nu, nu ) ) * lat( x - nu, nu );
	}

	Matrix u0 = lat( x, mu );
	Matrix u1 = mutate( u0, randGen );
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

template<class Matrix, class Mutator, class RandGen>
inline bool update0( LinkLattice<Matrix>& lat, double beta, Mutator mutate, RandGen& randGen ) {
	uniform_int<int> dist( 0, lat.size() - 1 );
	Site x(
		dist( randGen ),
		dist( randGen ),
		dist( randGen ),
		dist( randGen )
	);
	int mu = uniform_int<int>( 0, 3 )( randGen );

	return updateLocal( lat, x, mu, beta, mutate, randGen );
}

template<class Matrix, class Mutator, class RandGen>
inline void update1( LinkLattice<Matrix>& lat, double beta, Mutator mutate, RandGen& randGen ) {
	for( int x0 = 0; x0 < lat.size(); ++x0 )
	for( int x1 = 0; x1 < lat.size(); ++x1 )
	for( int x2 = 0; x2 < lat.size(); ++x2 )
	for( int x3 = 0; x3 < lat.size(); ++x3 ) {
		Site x( x0, x1, x2, x3 );
		for( int mu = 0; mu < 4; ++mu ) {
			return updateLocal( lat, x, mu, beta, mutate, randGen );
		}
	}
}