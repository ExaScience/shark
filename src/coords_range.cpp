/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <shark/coords_range.hpp>
#include <cassert>                    // assert

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim>
coords_range<ndim> coords_range<ndim>::overlap(coords_range<ndim> other) const {
	coords_range<ndim> r;
	seq<0,ndim>::for_each([this, &r, &other](int d) {
		r.lower[d] = std::max(lower[d], other.lower[d]);
		r.upper[d] = std::min(upper[d], other.upper[d]);
	});
	return r;
}

template<int ndim>
bool coords_range<ndim>::contains(coords<ndim> i) const {
	return seq<0,ndim>::all_of([this, &i](int d) {
		return lower[d] <= i[d] && upper[d] > i[d];
	});
}

template<int ndim>
bool coords_range<ndim>::contains(coords_range<ndim> other) const {
	return seq<0,ndim>::all_of([this, &other](int d) {
		return lower[d] <= other.lower[d] && upper[d] >= other.upper[d];
	});
}

template<int ndim>
coords<ndim+1> coords_range<ndim>::stride(coords<ndim> bw) const {
	coords<ndim+1> ld;
	ld[ndim] = 1;
	for(int d = ndim-1; d >= 0; d--)
		ld[d] = ld[d+1] * (upper[d] - lower[d] + 2 * bw[d]);
	return ld;
}


template<int ndim>
ostream& shark::ndim::operator<<(ostream& out, const coords_range<ndim>& r) {
	out << "[" << r.lower[0] << ":" << r.upper[0] << ")";
	seq<1,ndim-1>::for_each([&out,&r](int d) { out << " x " << "[" << r.lower[d] << ":" << r.upper[d] << ")"; });
	return out;
}

#if defined(SHARK_TBB_SCHED) && defined(SHARK_THREAD_BLOCK_DIST)

template<int ndim>
split_range<ndim>::split_range(coords_range<ndim> r, coord grainsize): r(r), grainsize(grainsize) {
}

template<int ndim>
bool split_range<ndim>::empty() const {
	return !seq<0,ndim>::all_of([this](int d) { return r.lower[d] < r.upper[d]; });
}

template<int ndim>
bool split_range<ndim>::is_divisible() const {
	return r.size() > grainsize;
}

template<int ndim>
void split_range<ndim>::split(split_range<ndim>& left, split_range<ndim>& right) const {
	assert(is_divisible());
	// Determine which dim to split
	coord size = coord();
	int ds = ndim;
	for(int d = ndim-1; d>=0; d--) {
		coord cand = r.upper[d] - r.lower[d];
		if(cand > size) {
			size = cand;
			ds = d;
		}
	}
	left.grainsize = grainsize;
	right.grainsize = grainsize;
	left.r = r;
	right.r = r;
	left.r.upper[ds] = r.lower[ds] + (r.upper[ds] - r.lower[ds]) / 2u;
	right.r.lower[ds] = left.r.upper[ds];
}

template<int ndim>
split_range<ndim>::split_range(split_range<ndim>& other, tbb::split): r(other.r), grainsize(other.grainsize) {
	split(*this, other);
}

#endif

// Set-up instantiations

#define SYMBD(d) template struct shark::ndim::coords_range<d>;
#include "inst_dim"
#undef SYMBD

#define SYMBD(d) template ostream& shark::ndim::operator<< <d>(ostream&, const coords_range<d>&);
#include "inst_dim"
#undef SYMBD

#if defined(SHARK_TBB_SCHED) && defined(SHARK_THREAD_BLOCK_DIST)
#define SYMBD(d) template struct shark::ndim::split_range<d>;
#include "inst_dim"
#undef SYMBD
#endif
