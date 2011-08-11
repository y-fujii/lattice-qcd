#pragma once

#include <iterator>
#include <vector>
#include <cmath>
#include <tr1/tuple>
#include "misc.hpp"
#include "link-lattice.hpp"

using namespace std;


struct GaugeAction {
	GaugeAction( double g ): _g( g ) {}

	template<class Lattice>
	double operator()( Lattice& lat, Site<Lattice::ndim> const& x, int mu, typename Lattice::Elem const& u ) {
		typename Lattice::Elem v = zero( typename Lattice::Elem() );
		for( int nu = 0; nu < Lattice::ndim; ++nu ) {
			if( mu == nu ) continue;
			Site<Lattice::ndim> y = x - nu;
			v = v + lat( x + mu, nu ) * inv( lat( x, nu ) * lat( x + nu, mu ) );
			v = v + inv( lat( y, mu ) * lat( y + mu, nu ) ) * lat( y, nu );
		}

		return (-_g) * real( ntrace( (u - lat( x, mu )) * v ) );
	}

	private:
		double const _g;
};

template<class Lattice, class ActionDiff, class Mutator, class RandGen>
double update( Lattice& lat, ActionDiff actionDiff, Mutator mutate, RandGen& rng ) {
	int nAcc = 0;
	Site<Lattice::ndim> x;
	x.assign( 0 );
	do {
		for( int mu = 0; mu < Lattice::ndim; ++mu ) {
			typename Lattice::Elem u = mutate( lat( x, mu ), rng );
			double dS = actionDiff( lat, x, mu, u );
			double p = uniform_real_ex<double>( 0.0, 1.0 )( rng );
			if( p <= exp( -dS ) ) {
				lat( x, mu ) = u;
				++nAcc;
			}
		}
	} while( lat.next( x ) );

	return double( nAcc ) / lat.nLinks();
}

template<class Lattice, class ActionDiff, class Mutator, class RandGen>
double update2( Lattice& lat, ActionDiff actionDiff, Mutator mutate, RandGen& rng ) {
	vector< tr1::tuple<Site<Lattice::ndim>, int> > xlinks;
	Site<Lattice::ndim> x;
	x.assign( 0 );
	do {
		for( int mu = 0; mu < Lattice::ndim; ++mu ) {
			xlinks.push_back( tr1::make_tuple( x, mu ) );
		}
	} while( lat.next( x ) );
	random_shuffle_ex( xlinks.begin(), xlinks.end(), rng );

	int nAcc = 0;
	for( size_t i = 0; i < xlinks.size(); ++i ) {
		Site<Lattice::ndim> x;
		int mu;
		tr1::tie( x, mu ) = xlinks[i];
		typename Lattice::Elem u = mutate( lat( x, mu ), rng );
		double dS = actionDiff( lat, x, mu, u );
		double p = uniform_real_ex<double>( 0.0, 1.0 )( rng );
		if( p <= exp( -dS ) ) {
			lat( x, mu ) = u;
			++nAcc;
		}
	}
	return double( nAcc ) / lat.nLinks();
}
