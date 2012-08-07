
#include <shark/domain.hpp>
#include <shark/expr.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim,typename T>
ConstExp<ndim,T>::ConstExp(const Domain<ndim>& dom, const T& val): dom(dom), val(val) { }

template<int ndim,typename T>
ConstExp<ndim,T> shark::ndim::constant(const Domain<ndim>& dom, const T& val) {
	return ConstExp<ndim,T>(dom, val);
}

template<int ndim,typename T>
ConstAcc<ndim,T>::ConstAcc(const ConstExp<ndim,T>& c): c(c) { }


// Set-up instantiations

#include "types"

#define SYMBH(d,T) template class ConstAcc<d,T>; 
#include "inst_dimtype_ho"
#undef SYMBH

#define SYMBH(d,T) template class ConstExp<d,T>; 
#include "inst_dimtype_ho"
#undef SYMBH

#define SYMBH(d,T) template ConstExp<d,T> shark::ndim::constant(const Domain<d>&, const T&); 
#include "inst_dimtype_ho"
#undef SYMBH
