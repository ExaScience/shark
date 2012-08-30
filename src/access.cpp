
#include <cassert>
#include <stdexcept>
#include <shark/access.hpp>
#include "comm_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim, typename T>
void Access<ndim,T>::reset() {
	ga = 0;
}

template<int ndim, typename T>
Access<ndim,T>::Access(): ga(0) { }

template<int ndim, typename T>
Access<ndim,T>::Access(Access<ndim,T>&& other): ga(other.ga), lower(other.lower) {
	other.reset();
}

template<int ndim,typename T>
Access<ndim,T>& Access<ndim,T>::operator=(Access<ndim,T>&& other) {
	if(*this)
		release();
	ga = other.ga;
	lower = other.lower;
	other.reset();
	return *this;
}

template<int ndim,typename T>
Access<ndim,T>::~Access() {
	if(*this)
		release();
}

template<int ndim,typename T>
Access<ndim,T>::Access(const GlobalArray<ndim,T>& ga): ga(&ga), lower(ga.domain().local().lower) {
	assert(ga);
#if defined(SHARK_MPI_COMM)
	if(ga.lc++ == 0)
		MPI_Win_lock(MPI_LOCK_SHARED, ga.domain().group.procid, 0, ga.impl->win);
#endif
}

template<int ndim,typename T>
Access<ndim,T>::Access(GlobalArray<ndim,T>& ga): ga(&ga), lower(ga.domain().local().lower) {
	assert(ga);
#if defined(SHARK_MPI_COMM)
	if(ga.lc++ == 0)
		MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ga.domain().group.procid, 0, ga.impl->win);
	else
		throw logic_error("Trying to create non-first mutable Access object");
#endif
}

template<int ndim,typename T>
void Access<ndim,T>::release() {
#if defined(SHARK_MPI_COMM)
	if(--ga->lc == 0)
		MPI_Win_unlock(ga->domain().group.procid, ga->impl->win);
#endif
}

// Set-up instantiations

#include "types"

#define SYMBDT(d,T) template class shark::ndim::Access<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

