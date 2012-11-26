/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <shark/coords.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim>
ostream& shark::ndim::operator<<(ostream& out, const coords<ndim>& i) {
	out << i[0];
	seq<1,ndim-1>::for_each([&out,&i](int d) { out << "," << i[d]; });
	return out;
}

template<int ndim>
array<size_t,ndim-1> shark::ndim::essential_lead(coords<ndim+1> ld) {
	array<size_t,ndim-1> eld;
	seq<1,ndim-1>::for_each([&eld,&ld](int d) { eld[d-1] = ld[d]; });
	return eld;
}

// Set-up instantiations

#define SYMBD(d) template ostream& shark::ndim::operator<< <d>(ostream&, const coords<d>&);
#include "inst_dim"
#undef SYMBD

#define SYMBD(d) template array<size_t,d-1> shark::ndim::essential_lead<d>(coords<d+1>);
#include "inst_dim"
#undef SYMBD
