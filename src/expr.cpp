/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <shark/domain.hpp>
#include <shark/expr.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim, typename T>
Const<ndim,T>::Const(const T& val): val(val) { }

template<int ndim, typename T>
Const<ndim,T>::~Const() { }

template<int ndim, typename T>
NullaryExp<ndim,Const<ndim,T>> shark::ndim::constant(const Domain<ndim>& dom, const T& val) {
	return NullaryExp<ndim,Const<ndim,T>>(dom, Const<ndim,T>(val));
}

template<int ndim, typename T>
NullaryExp<ndim,Const<ndim,T>> shark::ndim::constant(const Domain<ndim>& dom, coords_range<ndim> r, const T& val) {
	return NullaryExp<ndim,Const<ndim,T>>(dom, r, Const<ndim,T>(val));
}

template<int ndim, typename T>
Buffer<ndim,T>::Buffer(coords_range<ndim> r, coords<ndim+1> ld, const T* ptr): r(r), ld(ld), ptr(ptr) { }

template<int ndim, typename T>
Buffer<ndim,T>::~Buffer() { }

template<int ndim, typename T>
NullaryExp<ndim,Buffer<ndim,T>> shark::ndim::buffer(const Domain<ndim>& dom, coords_range<ndim> r, coords<ndim+1> ld, const T* ptr) {
	return NullaryExp<ndim,Buffer<ndim,T>>(dom, r, Buffer<ndim,T>(r, ld, ptr));
}

template<int ndim>
NullaryExp<ndim,CoordVec<ndim>> shark::ndim::coord_vec(const Domain<ndim>& dom) {
	return NullaryExp<ndim,CoordVec<ndim>>(dom, CoordVec<ndim>());
}

// Set-up instantiations

#include "types"

#define SYMBDT(d,T) template class shark::ndim::Const<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template NullaryExp<d,Const<d,T>> shark::ndim::constant(const Domain<d>&, const T&); 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template NullaryExp<d,Const<d,T>> shark::ndim::constant(const Domain<d>&, coords_range<d>, const T&); 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template class shark::ndim::Buffer<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template NullaryExp<d,Buffer<d,T>> shark::ndim::buffer(const Domain<d>&, coords_range<d>, coords<d+1>, const T*);
#include "inst_dimtype"
#undef SYMBDT

#define SYMBD(d) template NullaryExp<d,CoordVec<d>> shark::ndim::coord_vec(const Domain<d>&);
#include "inst_dim"
#undef SYMBD

