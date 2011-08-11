#pragma once

#include <complex>
#include <vector>
#include <tr1/random>
#include <Eigen/Core>
#include "misc.hpp"

using namespace std;
using namespace tr1;
using namespace Eigen;


template<class Derived>
inline typename Derived::PlainObject one( MatrixBase<Derived> const& ) {
	return Derived::Identity();
}

template<class Derived>
inline typename Derived::PlainObject zero( MatrixBase<Derived> const& ) {
	return Derived::Zero();
}

template<class Derived>
inline complex<double> ntrace( MatrixBase<Derived> const& m ) {
	return m.trace() * (1.0 / m.rows());
}

template<class Derived>
inline typename Derived::PlainObject inv( MatrixBase<Derived> const& m ) {
	return m.adjoint(); // * (N / m.squaredNorm());
}

template<class Derived>
inline typename Derived::PlainObject orthonormalize( MatrixBase<Derived> const& m ) {
	typename Derived::PlainObject n;
	for( size_t i = 0; i < m.rows(); ++i ) {
		typename internal::plain_row_type<Derived>::type v = m.row( i );
		for( size_t j = 0; j < i; ++j ) {
			//v -= n.row( j ) * n.row( j ).dot( n.row( i ) );
			v -= n.row( j ) * n.row( j ).dot( v );
		}
		n.row( i ) = v.normalized();
	}
	return n;
}

template<class Derived, class RandGen>
inline typename Derived::PlainObject orthonormalize( MatrixBase<Derived> const& m, RandGen& rng ) {
	vector<size_t> is( m.rows() );
	for( size_t i = 0; i < is.size(); ++i ) {
		is[i] = i;
	}
	random_shuffle_ex( is.begin(), is.end(), rng );

	typename Derived::PlainObject n;
	for( size_t i = 0; i < m.rows(); ++i ) {
		typename internal::plain_row_type<Derived>::type v = m.row( is[i] );
		for( size_t j = 0; j < i; ++j ) {
			v -= n.row( is[j] ) * n.row( is[j] ).dot( v );
		}
		n.row( is[i] ) = v.normalized();
	}
	return n;
}

struct Greater1st {
	template<class T, class U>
	bool operator()( pair<T, U> const& x, pair<T, U> const& y ) {
		return x.first > y.first;
	}
};

template<class Derived>
inline typename Derived::PlainObject orthonormalize_ex( MatrixBase<Derived> const& m ) {
	vector< pair<double, size_t> > ix( m.rows() );
	for( size_t i = 0; i < ix.size(); ++i ) {
		ix[i] = make_pair( m.row( i ).squaredNorm(), i );
	}
	sort( ix.begin(), ix.end(), Greater1st() );

	typename Derived::PlainObject n;
	for( size_t i = 0; i < m.rows(); ++i ) {
		size_t mi = ix[i].second;
		typename internal::plain_row_type<Derived>::type v = m.row( mi );
		for( size_t j = 0; j < i; ++j ) {
			size_t mj = ix[j].second;
			v -= n.row( mj ) * n.row( mj ).dot( v );
		}
		n.row( mi ) = v.normalized();
	}
	return n;
}
