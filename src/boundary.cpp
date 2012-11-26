/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <shark/boundary.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim, typename T>
Boundary<ndim,T>::type::~type() { }

template<int ndim, typename T>
typename Boundary<ndim,T>::type* Boundary<ndim,T>::type::clone() const {
	return new type();
}

template<int ndim, typename T>
Boundary<ndim,T>::periodic_type::periodic_type() { }

template<int ndim, typename T>
Boundary<ndim,T>::periodic_type::~periodic_type() { }

template<int ndim, typename T>
typename Boundary<ndim,T>::periodic_type* Boundary<ndim,T>::periodic_type::clone() const {
	return new periodic_type();
}

template<int ndim, typename T>
Boundary<ndim,T>::Boundary(type* t): t(t) { }

template<int ndim, typename T>
Boundary<ndim,T>::Boundary(): t(new type()) { }

template<int ndim, typename T>
array<Boundary<ndim,T>,ndim> Boundary<ndim,T>::all_unmanaged() {
	return array<Boundary<ndim,T>,ndim>();
}

template<int ndim, typename T>
Boundary<ndim,T>::~Boundary() { }

template<int ndim, typename T>
Boundary<ndim,T>::Boundary(const Boundary<ndim,T>& other): t(other.t->clone()) { }

template<int ndim, typename T>
Boundary<ndim,T>& Boundary<ndim,T>::operator=(const Boundary<ndim,T>& other) {
	t.reset(other.t->clone());
	return *this;
}

template<int ndim, typename T>
Boundary<ndim,T> Boundary<ndim,T>::periodic() {
	return Boundary<ndim,T>(new periodic_type());
}

template<int ndim, typename T>
array<Boundary<ndim,T>,ndim> Boundary<ndim,T>::all_periodic() {
	array<Boundary<ndim,T>,ndim> bds;
	for(int d = 0; d < ndim; d++)
		bds[d].t.reset(new periodic_type());
	return bds;
}

template<int ndim, typename T>
Boundary<ndim,T> Boundary<ndim,T>::constant(const T& val) {
	return fixed([val](coords<ndim>) {
		return val;
	});
}

template<int ndim, typename T>
array<Boundary<ndim,T>,ndim> Boundary<ndim,T>::all_constant(const T& val) {
	return all_fixed([val](coords<ndim>) {
		return val;
	});
}

// Set-up instantiations

#include "types"

#define SYMBDT(d,T) template class shark::ndim::Boundary<d,T>; 
#include "inst_dimtype"
#undef SYMBDT
