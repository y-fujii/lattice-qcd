#pragma once

#include <vector>
#include <cassert>


struct Site {
	Site() {}

	Site( int x0, int x1, int x2, int x3 ) {
		_arr[0] = x0;
		_arr[1] = x1;
		_arr[2] = x2;
		_arr[3] = x3;
	}

	int& operator[]( int i ) {
		return _arr[i];
	}

	int operator[]( int i ) const {
		return _arr[i];
	}

	private:
		int _arr[4];
};

inline Site operator+( Site const& x, int mu ) {
	assert( 0 <= mu && mu < 4 );
	Site z = x;
	z[mu] += 1;
	return z;
}

inline Site operator-( Site const& x, int mu ) {
	assert( 0 <= mu && mu < 4 );
	Site z = x;
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

template<class T>
struct LinkLattice {
	LinkLattice( int n ):
		_size( n ),
		_array( n * n * n * n * 4 )
	{}

	T& operator()( Site const& x, int mu ) {
		assert( 0 <= mu && mu < 4 );
		int idx = 0;
		for( int i = 0; i < 4; ++i ) {
			idx = idx * _size + periodic( x[i], _size );
		}
		return _array[idx * 4 + mu];
	}

	T& operator()( int i, int mu ) {
		assert( 0 <= mu && mu < 4 );
		assert( 0 <= i && i < nSites() );
		return _array[i * 4 + mu];
	}

	T& operator()( int i ) {
		assert( 0 <= i && i < nLinks() );
		return _array[i];
	}

	int nSites() const {
		return _size * _size * _size * _size;
	}

	int nLinks() const {
		return nSites() * 4;
	}

	int size() const {
		return _size;
	}

	private:
		int const _size;
		std::vector<T> _array;
};
