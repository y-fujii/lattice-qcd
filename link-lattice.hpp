#pragma once

#include <vector>
#include <cassert>


template<int D = 4>
struct Site {
	static inline Site zero() {
		Site<D> x;
		for( int i = 0; i < D; ++i ) {
			x[i] = 0;
		}
		return x;
	}

	Site() {}

	int& operator[]( int i ) {
		assert( 0 <= i && i < D );
		return _arr[i];
	}

	int operator[]( int i ) const {
		assert( 0 <= i && i < D );
		return _arr[i];
	}

	private:
		int _arr[D];
};

template<int D>
inline Site<D> operator+( Site<D> const& x, int mu ) {
	assert( 0 <= mu && mu < D );
	Site<D> z = x;
	z[mu] += 1;
	return z;
}

template<int D>
inline Site<D> operator-( Site<D> const& x, int mu ) {
	assert( 0 <= mu && mu < D );
	Site<D> z = x;
	z[mu] -= 1;
	return z;
}

inline int periodic( int x, int N ) {
	assert( 0 <= x + N && x - N < N );
	if( x < 0 ) {
		return x + N;
	}
	else if( x < N ) {
		return x;
	}
	else {
		return x - N;
	}
}

template<int n>
inline int constPow( int x ) {
	return x * constPow<n - 1>( x );
}

template<>
inline int constPow<0>( int x ) {
	return 1;
}

template<class T, int D = 4>
struct LinkLattice {
	LinkLattice( int n ):
		_size( n ),
		_array( constPow<D>( n ) * D )
	{}

	T& operator()( Site<D> const& x, int mu ) {
		assert( 0 <= mu && mu < D );
		int idx = 0;
		for( int i = 0; i < D; ++i ) {
			idx = idx * _size + periodic( x[i], _size );
		}
		return _array[idx * D + mu];
	}

	T& operator()( int i, int mu ) {
		assert( 0 <= mu && mu < D );
		assert( 0 <= i && i < nSites() );
		return _array[i * D + mu];
	}

	T& operator()( int i ) {
		assert( 0 <= i && i < nLinks() );
		return _array[i];
	}

	int nSites() const {
		return constPow<D>( _size );
	}

	int nLinks() const {
		return nSites() * D;
	}

	int size() const {
		return _size;
	}

	bool next( Site<D>& x ) {
		for( int i = 0; i < D; ++i ) {
			if( ++x[i] < _size ) {
				return true;
			}
			x[i] = 0;
		}
		return false;
	}

	private:
		int const _size;
		std::vector<T> _array;
};
