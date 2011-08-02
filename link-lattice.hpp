#pragma once

#include <tr1/array>
#include <vector>
#include <cassert>


template<int D = 4>
struct Site: std::tr1::array<int, D> {
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

template<int n>
inline int constPow( int x ) {
	return x * constPow<n - 1>( x );
}

template<>
inline int constPow<0>( int ) {
	return 1;
}

template<class T, int D = 4>
struct LinkLattice {
	typedef typename std::vector<T>::iterator iterator;
	typedef T Elem;
	static int const ndim = D;

	LinkLattice( int n ):
		_size( n ),
		_array( constPow<D>( n ) * D )
	{}

	T& operator()( Site<D> const& x, int mu ) {
		assert( 0 <= mu && mu < D );
		int idx = 0;
		for( int i = 0; i < D; ++i ) {
			assert( 0 <= x[i] + _size && x[i] - _size < _size );
			idx = idx * _size + x[i];
			if( __builtin_expect( x[i] < 0, false ) ) {
				idx += _size;
			}
			else if( __builtin_expect( _size <= x[i], false ) ) {
				idx -= _size;
			}
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

	int len() const {
		return _size;
	}

	bool next( Site<D>& x ) const {
		for( int i = D - 1; i >= 0; --i ) {
			if( __builtin_expect( ++x[i] < _size, true ) ) {
				return true;
			}
			x[i] = 0;
		}
		return false;
	}

	iterator begin() {
		return _array.begin();
	}

	iterator end() {
		return _array.end();
	}

	private:
		int const _size;
		std::vector<T> _array;
};
