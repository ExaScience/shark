
#include <cassert>                     // assert

#include <shark/domain.hpp>
#include <shark/globalarray.hpp>
#include "mpi_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim, typename T>
const int GlobalArray<ndim,T>::number_of_dimensions;

template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(): ptr(0) { }

template<int ndim, typename T>
GlobalArray<ndim,T>::~GlobalArray() {
	if(*this) {
		assert(lc == 0);
		deallocate();
	}
}

template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(GlobalArray<ndim,T>&& other):
	dom(other.dom),
	gw(other.gw),
	gc(other.gc),
	ptr(other.ptr),
	impl(std::move(other.impl)),
	ld(std::move(other.ld)),
	lc(0)
{
	assert(!other || other.lc == 0);
	other.reset();
}

template<int ndim, typename T>
GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator=(GlobalArray<ndim,T>&& other) {
	assert(!other || other.lc == 0);
	if(*this)
		deallocate();
	dom = other.dom;
	gw = other.gw;
	gc = other.gc;
	ptr = other.ptr;
	impl = std::move(other.impl);
	ld = std::move(other.ld);
	lc = 0;
	other.reset();
	return *this;
}

template <int ndim, typename T>
void GlobalArray<ndim,T>::reset() {
	ptr = 0;
}

template<int ndim,typename T>
coords<ndim+1> GlobalArray<ndim,T>::stride(coords<ndim> count, coords<ndim> ghost_width) {
	coords<ndim+1> ld;
	ld[ndim] = 1;
	for(int d = ndim-1; d >= 0; d--)
		ld[d] = ld[d+1] * (count[d] + 2 * ghost_width[d]);
	return ld;
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::allocate(const Domain<ndim>& domain, coords<ndim> ghost_width, bool ghost_corners, T **ptr, GAImpl<ndim>& impl, coords<ndim+1>& ld) {
	coords<ndim> count = domain.count();

	ld = stride(count, ghost_width);

	// Create ghost types
	MPI_Aint lb, extent;
	MPI_Type_get_extent(mpi_type<T>::t, &lb, &extent);
	for(int di = 0; di < ndim; di++) {
		MPI_Type_dup(mpi_type<T>::t, &impl.ghost[di]);
		for(int d = ndim-1; d >= 0; d--) {
			MPI_Datatype tmp = impl.ghost[di];
			coord n =
				d == di ?
				ghost_width[d] :
				(ghost_corners && d < di ?
					count[d] + 2 * ghost_width[d] :
					count[d]);
			MPI_Type_create_hvector(n, 1, ld[d+1]*extent, tmp, &impl.ghost[di]);
			MPI_Type_free(&tmp);
		}
		MPI_Type_commit(&impl.ghost[di]);
	}

	// Allocate memory
	MPI_Aint size = ld[0] * sizeof(T);
	MPI_Alloc_mem(size, MPI_INFO_NULL, ptr);

	// Create window
	MPI_Win_create(*ptr, size, static_cast<int>(sizeof(T)), MPI_INFO_NULL, domain.group.impl->comm, &impl.win);
}

// Constructors
template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(const Domain<ndim>& domain, coords<ndim> ghost_width, bool ghost_corners):
  dom(&domain), gw(ghost_width), gc(ghost_corners), impl(new GAImpl<ndim>()), lc(0) {
	allocate(domain, ghost_width, ghost_corners, &ptr, *impl, ld);
}

// Destructor

template<int ndim,typename T>
void GlobalArray<ndim,T>::deallocate() {
	for(int di = 0; di < ndim; di++)
		MPI_Type_free(&impl->ghost[di]);

	MPI_Win_free(&impl->win);
	MPI_Free_mem(ptr);
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::reshape(const Domain<ndim>& domain) {
	assert(this->domain().equiv(domain));

	T* ptr;
	unique_ptr<GAImpl<ndim>> impl(new GAImpl<ndim>());
	coords<ndim+1> ld;

	allocate(domain, ghost_width(), ghost_corners(), &ptr, *impl, ld);
	if(domain.hasData()) {
		array<size_t,ndim-1> eld;
		size_t off = gw[ndim-1];
		for(int d = ndim-2; d >= 0; d--) {
			eld[d] = ld[d+1];
			off += gw[d] * eld[d];
		}
		this->get(domain.local(), eld, ptr+off);
	}

	deallocate();
	this->dom = &domain;
	this->ptr = ptr;
	this->impl = move(impl);
	this->ld = move(ld);
}

// Inspectors

template<int ndim,typename T>
void GlobalArray<ndim,T>::update() const {
	bool any = false;
	const coords<ndim> gw = ghost_width();
	const bool gc = ghost_corners();
	const coords<ndim> count = domain().count();

	for(int di = 0; di < ndim; di++) {
		if(gw[di] > 0) {
			any = true;
			const int prev = domain().shiftd(di, -1);
			const int next = domain().shiftd(di,  1);
			coords<ndim> fronti, backi;
			// backward
			for(int d = 0; d < ndim; d++) {
				fronti[d] = gc && d < di ? -gw[d] : 0;
				backi[d] = d == di ? count[d] : fronti[d];
			}
			MPI_Sendrecv(&access(fronti), 1, impl->ghost[di], prev, 2*di,
					&access(backi), 1, impl->ghost[di], next, 2*di,
					domain().group.impl->comm, MPI_STATUS_IGNORE);
			// forward
			fronti[di] -= gw[di];
			backi[di] -= gw[di];
			MPI_Sendrecv(&access(backi), 1, impl->ghost[di], next, 2*di+1,
					&access(fronti), 1, impl->ghost[di], prev, 2*di+1,
					domain().group.impl->comm, MPI_STATUS_IGNORE);
		}
	}

	if(!any)
		domain().sync();
	/* 
		int recvcounts[nprocs];
		int displs[nprocs];
		for(int k = 0; k < nprocs; k++) {
			recvcounts[k] = row_domain().count(k) * cols();
			displs[k] = row_domain().lower(k) * cols();
		}
		MPI_Allgatherv(MPI_IN_PLACE, 0, mpi_type<T>(), ptr, recvcounts, displs, mpi_type<T>(), MPI_COMM_WORLD);
	*/
}

template<int ndim,typename T>
class GlobalArray<ndim,T>::RMAOp {
	typedef array<int,ndim> pcoords;

	const Domain<ndim>& dom;
	const coords_range<ndim> range;
	const coords<ndim> gw;
	const array<size_t,ndim-1> ognld;
	pcoords lip, uip;
	coords<ndim> loff, uoff;

	template<int d,typename Op>
	INLINE typename enable_if<d < ndim>::type opd(const Op& op, pcoords& ip, coords_range<ndim>& i);

	template<int d,typename Op>
	INLINE typename enable_if<d == ndim>::type opd(const Op& op, pcoords& ip, coords_range<ndim>& i);

public:
	RMAOp(const Domain<ndim>& dom, coords_range<ndim> range, coords<ndim> gw, array<size_t,ndim-1> ognld);

	template<typename Op>
	void op(Op op);
};

template<int ndim,typename T> template<int d,typename Op>
inline typename enable_if<d < ndim>::type GlobalArray<ndim,T>::RMAOp::opd(const Op& op, pcoords& ip, coords_range<ndim>& i) {
	for(ip[d] = lip[d]; ip[d] <= uip[d]; ip[d]++) {
		i.lower[d] = ip[d] == lip[d] ? loff[d] : 0;
		i.upper[d] = ip[d] == uip[d] ? uoff[d] : dom.nd[d][ip[d]+1] - dom.nd[d][ip[d]];
		opd<d+1>(op, ip, i);
	}
}

template<int ndim,typename T> template<int d,typename Op>
inline typename enable_if<d == ndim>::type GlobalArray<ndim,T>::RMAOp::opd(const Op& op, pcoords& ip, coords_range<ndim>& i) {
	// Target
	int id = dom.pindex(ip);
	coords<ndim+1> tgtld = stride(dom.count(id), gw);

	// Offsets
	size_t ognoff = 0;
	MPI_Aint tgtoff = 0;
	for(int di = 0; di < ndim; di++) {
		ognoff += (di == ndim-1 ? 1 : ognld[di]) * (dom.nd[di][ip[di]] + i.lower[di] - range.lower[di]);
		tgtoff += tgtld[di+1] * (gw[di] + i.lower[di]);
	}

	// Types
	MPI_Datatype ognt, tgtt;
	{
		MPI_Aint lb, extent;
		MPI_Type_get_extent(mpi_type<T>::t, &lb, &extent);

		MPI_Type_dup(mpi_type<T>::t, &ognt);
		MPI_Type_dup(mpi_type<T>::t, &tgtt);
		for(int di = ndim-1; di >= 0; di--) {
			MPI_Datatype tmp;
			coord n = i.upper[di] - i.lower[di];
			tmp = ognt;
			MPI_Type_create_hvector(n, 1, di == ndim-1 ? extent : ognld[di]*extent, tmp, &ognt);
			MPI_Type_free(&tmp);
			tmp = tgtt;
			MPI_Type_create_hvector(n, 1, tgtld[di+1]*extent, tmp, &tgtt);
			MPI_Type_free(&tmp);
		}
		MPI_Type_commit(&ognt);
		MPI_Type_commit(&tgtt);
	}
	op(id, ognoff, ognt, tgtoff, tgtt);

	MPI_Type_free(&tgtt);
	MPI_Type_free(&ognt);
}

template<int ndim,typename T>
GlobalArray<ndim,T>::RMAOp::RMAOp(const Domain<ndim>& dom, coords_range<ndim> range, coords<ndim> gw, array<size_t,ndim-1> ognld)
: dom(dom), range(range), gw(gw), ognld(ognld) {
	assert(range.lower >= coords<ndim>() && range.lower < dom.n);
	assert(range.upper >= coords<ndim>() && range.upper <= dom.n);
	assert(range.lower <= range.upper);
	dom.find(range.lower, lip, loff);
	dom.find(range.upper, uip, uoff);
	for(int d = 0; d < ndim; d++)
		if(uoff[d] == 0) {
			uip[d]--;
			uoff[d] = dom.nd[d][uip[d]+1] - dom.nd[d][uip[d]];
		}
}

template<int ndim,typename T> template<typename Op>
void GlobalArray<ndim,T>::RMAOp::op(Op op) {
	pcoords ip;
	coords_range<ndim> i;
	opd<0>(op, ip, i);
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::get(coords_range<ndim> range, array<size_t,ndim-1> ld, T* buf) const {
	RMAOp(domain(), range, ghost_width(), ld).op(
		[&impl,buf](int id, size_t ognoff, MPI_Datatype ognt, MPI_Aint tgtoff, MPI_Datatype tgtt) {
			MPI_Win_lock(MPI_LOCK_SHARED, id, 0, impl->win);
			MPI_Get(buf+ognoff, 1, ognt, id, tgtoff, 1, tgtt, impl->win);
			MPI_Win_unlock(id, impl->win);
		}
	);
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::put(coords_range<ndim> range, array<size_t,ndim-1> ld, const T* buf) {
	RMAOp(domain(), range, ghost_width(), ld).op(
		[&impl,buf](int id, size_t ognoff, MPI_Datatype ognt, MPI_Aint tgtoff, MPI_Datatype tgtt) {
			MPI_Win_lock(MPI_LOCK_EXCLUSIVE, id, 0, impl->win);
			MPI_Put(const_cast<T*>(buf+ognoff), 1, ognt, id, tgtoff, 1, tgtt, impl->win);
			MPI_Win_unlock(id, impl->win);
		}
	);
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::accumulate(coords_range<ndim> range, array<size_t,ndim-1> ld, const T* buf) {
	RMAOp(domain(), range, ghost_width(), ld).op(
		[&impl,buf](int id, size_t ognoff, MPI_Datatype ognt, MPI_Aint tgtoff, MPI_Datatype tgtt) {
			MPI_Win_lock(MPI_LOCK_SHARED, id, 0, impl->win);
			MPI_Accumulate(const_cast<T*>(buf+ognoff), 1, ognt, id, tgtoff, 1, tgtt, MPI_SUM, impl->win);
			MPI_Win_unlock(id, impl->win);
		}
	);
}

// Set-up instantiations

#define SYMB(d,T) template class GlobalArray<d,T>; 
#include "inst_dimtype"
#undef SYMB
