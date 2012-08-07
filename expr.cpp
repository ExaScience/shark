
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

// Set-up instantiations

#include "types"

#define SYMBH(d,T) template class Const<d,T>; 
#include "inst_dimtype_ho"
#undef SYMBH

#define SYMBH(d,T) template NullaryExp<d,Const<d,T>> shark::ndim::constant(const Domain<d>&, const T&); 
#include "inst_dimtype_ho"
#undef SYMBH
