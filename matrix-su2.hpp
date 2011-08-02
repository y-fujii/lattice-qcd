#pragma once

#include <limits>


template<class T>
struct MatrixSU2 {
	static inline MatrixSU2 zero() {
		return MatrixSU2( 0.0, 0.0, 0.0, 0.0 );
	}

	static inline MatrixSU2 one() {
		return MatrixSU2( 1.0, 0.0, 0.0, 0.0 );
	}

#if defined( NDEBUG )
	MatrixSU2() {}
#else
	MatrixSU2():
		a0( std::numeric_limits<T>::signaling_NaN() ),
		a1( std::numeric_limits<T>::signaling_NaN() ),
		a2( std::numeric_limits<T>::signaling_NaN() ),
		a3( std::numeric_limits<T>::signaling_NaN() )
	{}
#endif

	MatrixSU2( T const& a0_, T const& a1_, T const& a2_, T const& a3_ ):
		a0( a0_ ), a1( a1_ ), a2( a2_ ), a3( a3_ ) {}

	T a0, a1, a2, a3;
};

template<class T>
inline T det( MatrixSU2<T> const& x ) {
	return x.a0 * x.a0 + x.a1 * x.a1 + x.a2 * x.a2 + x.a3 * x.a3;
}

template<class T>
inline T tr( MatrixSU2<T> const& x ) {
	return x.a0 + x.a0;
}

template<class T>
inline MatrixSU2<T> inv( MatrixSU2<T> const& x ) {
	return MatrixSU2<T>( x.a0, -x.a1, -x.a2, -x.a3 );
}

template<class T>
inline MatrixSU2<T> operator-( MatrixSU2<T> const& x ) {
	return MatrixSU2<T>( -x.a0, -x.a1, -x.a2, -x.a3 );
}

template<class T>
inline MatrixSU2<T> operator+( MatrixSU2<T> const& x, MatrixSU2<T> const& y ) {
	return MatrixSU2<T>( x.a0 + y.a0, x.a1 + y.a1, x.a2 + y.a2, x.a3 + y.a3 );
}

template<class T>
inline MatrixSU2<T> operator-( MatrixSU2<T> const& x, MatrixSU2<T> const& y ) {
	return MatrixSU2<T>( x.a0 - y.a0, x.a1 - y.a1, x.a2 - y.a2, x.a3 - y.a3 );
}

template<class T, class U>
inline MatrixSU2<T> operator*( MatrixSU2<T> const& x, U const& k ) {
	return MatrixSU2<T>( x.a0 * k, x.a1 * k, x.a2 * k, x.a3 * k );
}

template<class T>
inline MatrixSU2<T> operator*( MatrixSU2<T> const& x, MatrixSU2<T> const& y ) {
	return MatrixSU2<T>(
		x.a0 * y.a0 - x.a1 * y.a1 - x.a2 * y.a2 - x.a3 * y.a3,
		x.a1 * y.a0 + x.a0 * y.a1 - x.a2 * y.a3 + x.a3 * y.a2,
		x.a2 * y.a0 + x.a0 * y.a2 - x.a3 * y.a1 + x.a1 * y.a3,
		x.a3 * y.a0 + x.a0 * y.a3 - x.a1 * y.a2 + x.a2 * y.a1
	);
}
