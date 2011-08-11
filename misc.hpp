#pragma once

#include <cassert>
#include <cmath>


template<class T>
struct uniform_real_ex {
	uniform_real_ex( T const& min, T const& max ):
		_min( min ), _max( max ) {}

	template<class RandGen>
	T operator()( RandGen& gen ) {
		assert( gen.min() == 0 );
		return gen() * ((_max - _min) / (1.0 + gen.max())) + _min;
	}

	private:
		T const _min, _max;
};

template<class T>
struct normal_distribution_ex {
	normal_distribution_ex( T const& m, T const& s ):
		_mean( m ), _sigma( s ) {}

	template<class RandGen>
	double operator()( RandGen& rng ) {
		assert( rng.min() == 0 );

		T u = (rng() + 1) * (1.0 / (rng.max() + 1));
		T v = rng() * (2.0 * M_PI / (rng.max() + 1));
		T p = std::sqrt( -2.0 * std::log( u ) ) * std::sin( v );

		return p * _sigma + _mean;
	}

	private:
		T const _mean;
		T const _sigma;
};

template<class Iter, class RandGen>
inline void random_shuffle_ex( Iter bgn, Iter end, RandGen& rng ) {
	for( Iter it = bgn; it < end; ++it ) {
		// it <= jt < end
		Iter jt = it + std::tr1::uniform_int<size_t>( 0, (end - 1) - it )( rng );
		std::swap( *it, *jt );
	}
}

inline double real( double x ) {
		return x;
}

template<int n>
inline int constPow( int x ) {
	return x * constPow<n - 1>( x );
}

template<>
inline int constPow<0>( int ) {
	return 1;
}
