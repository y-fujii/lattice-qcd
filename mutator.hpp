#pragma once

#include <cmath>
#include <tr1/random>
#include "misc.hpp"

using namespace std;
using namespace tr1;


struct U1Mutator {
	U1Mutator( double e ): _eps( e ) {}

	template<class RandGen>
	complex<double> operator()( complex<double> const& u, RandGen& rng ) {
		double im = uniform_real_ex<double>( -_eps, +_eps )( rng );
		return complex<double>( cos( im ), sin( im ) ) * u;
	}

	private:
		double const _eps;
};

template<class RandGen>
inline MatrixSU2<double> randSU2( double eps, RandGen& rng ) {
	double x, y, z, r2;
	do {
		x = uniform_real_ex<double>( -eps, +eps )( rng );
		y = uniform_real_ex<double>( -eps, +eps )( rng );
		z = uniform_real_ex<double>( -eps, +eps )( rng );
		r2 = x * x + y * y + z * z;
	} while( r2 > eps * eps );

	return MatrixSU2<double>( sqrt( 1.0 - r2 ), x, y, z );
}

struct SU2Mutator {
	SU2Mutator( double e ): _eps( e ) {}

	template<class RandGen>
	MatrixSU2<double> operator()( MatrixSU2<double> const& u, RandGen& rng ) {
		return randSU2( _eps, rng ) * u;
	}

	private:
		double const _eps;
};

template<class Matrix, class RandGen>
inline Matrix randAntiHermite( double eps, RandGen& rng ) {
	Matrix m;
	double tr = 0.0;
	for( int i = 0; i < m.rows(); ++i ) {
		for( int j = 0; j < i; ++j ) {
			double re = normal_distribution_ex<double>( 0.0, eps )( rng );
			double im = normal_distribution_ex<double>( 0.0, eps )( rng );
			m( i, j ) = complex<double>( +re, im );
			m( j, i ) = complex<double>( -re, im );
		}
		double im = normal_distribution_ex<double>( 0.0, eps )( rng );
		m( i, i ) = complex<double>( 0.0, im );
		tr += im;
	}
	for( int i = 0; i < m.rows(); ++i ) {
		m( i, i ).imag() -= tr / m.rows();
	}
	return m;
}

template<class Matrix, class RandGen>
inline Matrix randSUn( double eps, RandGen& rng ) {
	Matrix m = Matrix::Identity();
	assert( m.rows() == m.cols() );
	int i, j;
	do {
		i = uniform_int<int>( 0, m.rows() - 1 )( rng );
		j = uniform_int<int>( 0, m.cols() - 1 )( rng );
	} while( i == j );

	MatrixSU2<double> su2 = randSU2( eps, rng );

	m( i, i ) = complex<double>( +su2.a0, +su2.a3 );
	m( i, j ) = complex<double>( -su2.a2, +su2.a1 );
	m( j, i ) = complex<double>( +su2.a2, +su2.a1 );
	m( j, j ) = complex<double>( +su2.a0, -su2.a3 );

	return m;
}

struct SUnMutator {
	SUnMutator( double e ): _eps( e ) {}

	template<class Derived, class RandGen>
	typename Derived::PlainObject operator()( MatrixBase<Derived> const& u, RandGen& rng ) {
		typedef typename Derived::PlainObject Matrix;
		/*
		Matrix m = orthonormalize_ex(
			Derived::Identity() + randAntiHermite<Matrix>( _eps, rng )
		);
		return m * u;
		*/
		/*
		Matrix a = randAntiHermite<Matrix>( _eps, rng );
		Matrix m = a * (
			(1.0 / 2.0) * a * (
				(1.0 / 3.0) * a * (
					(1.0 / 4.0) * a + Derived::Identity()
				).eval() + Derived::Identity()
			).eval() + Derived::Identity()
		).eval() + Derived::Identity();

		return orthonormalize_ex( m ) * u;
		*/
		return randSUn<Matrix>( _eps, rng ) * u;
	}

	private:
		double const _eps;
};
