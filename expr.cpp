
#include <shark/domain.hpp>
#include <shark/expr.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim, typename T>
Const<ndim,T>::Const(const T& val): val(val) { }

template<int ndim, typename T>
Const<ndim,T>::~Const() { }

template<int ndim,typename T>
NullaryExp<ndim,Const<ndim,T>> shark::ndim::constant(const Domain<ndim>& dom, const T& val) {
	return NullaryExp<ndim,Const<ndim,T>>(dom, Const<ndim,T>(val));
}

template<int ndim,typename T>
NullaryExp<ndim,Const<ndim,T>> shark::ndim::constant(const Domain<ndim>& dom, coords_range<ndim> r, const T& val) {
	return NullaryExp<ndim,Const<ndim,T>>(dom, r, Const<ndim,T>(val));
}

template<int ndim>
Coord<ndim>::Coord() { }

template<int ndim>
Coord<ndim>::~Coord() { }

template<int ndim>
NullaryExp<ndim,Coord<ndim>> shark::ndim::coord_vec(const Domain<ndim>& dom) {
	return NullaryExp<ndim,Coord<ndim>>(dom, Coord<ndim>());
}

// Set-up instantiations

#include "types"

#define SYMBDT(d,T) template class Const<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template NullaryExp<d,Const<d,T>> shark::ndim::constant(const Domain<d>&, const T&); 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template NullaryExp<d,Const<d,T>> shark::ndim::constant(const Domain<d>&, coords_range<d>, const T&); 
#include "inst_dimtype"
#undef SYMBDT

#define SYMBD(d) template class Coord<d>; 
#include "inst_dim"
#undef SYMBD

#define SYMBD(d) template NullaryExp<d,Coord<d>> shark::ndim::coord_vec(const Domain<d>&);
#include "inst_dim"
#undef SYMBD
