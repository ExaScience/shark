/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <iostream>
#include <shark/access.hpp>
#include <shark/boundary.hpp>
#include <shark/domain.hpp>
#include <shark/globalarray.hpp>
#include <shark/sparsearray.hpp>
#include "comm_impl.hpp"

#include <utility>                     // std::pair

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim, typename T>
const int GlobalArray<ndim,T>::number_of_dimensions;

template<int ndim, typename T>
ostream& GlobalArray<ndim,T>::log_out() const
{
	return *shark::log_out << domain().group.procid << ": " << this << ": ";
}

template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(): dom(0) { }

template<int ndim, typename T>
GlobalArray<ndim,T>::~GlobalArray()
{
	if(*this)
	{
		assert(lc == 0);
		deallocate();
	}
}

template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(GlobalArray<ndim,T>&& other):
	dom(other.dom),
	gw(other.gw),
	gc(other.gc),
	bd(std::move(other.bd)),
	ptr(other.ptr),
	impl(std::move(other.impl)),
	ld(std::move(other.ld)),
	ghost_back(std::move(other.ghost_back)),
	ghost_front(std::move(other.ghost_front)),
	lc(0)
{
	assert(!other || other.lc == 0);
	other.reset();
}

template<int ndim, typename T>
GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator=(GlobalArray<ndim,T>&& other)
{
	assert(!other || other.lc == 0);
	if(*this)
		deallocate();
	dom = other.dom;
	gw = other.gw;
	gc = other.gc;
	bd = std::move(other.bd);
	ptr = other.ptr;
	impl = std::move(other.impl);
	ld = std::move(other.ld);
	ghost_back = std::move(other.ghost_back);
	ghost_front = std::move(other.ghost_front);
	lc = 0;
	other.reset();
	return *this;
}

template<int ndim, typename T>
GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator=(const GlobalArray<ndim,T>& other) {
	region(domain().total()) = other;
	return *this;
}

template <int ndim, typename T>
void GlobalArray<ndim,T>::reset() {
	dom = 0;
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::allocate()
{
	const coords_range<ndim> local = domain().local();
	const coords<ndim> gw = ghost_width();
	const bool gc = ghost_corners(); 

	ld = local.stride(gw);

	for(int di = 0; di < ndim; di++)
		for(int d = 0; d < ndim; d++)
			if(d == di)
			{
				ghost_back [di].lower[d] = local.lower[d] - gw[d];
				ghost_back [di].upper[d] = local.lower[d];
				ghost_front[di].lower[d] = local.upper[d];
				ghost_front[di].upper[d] = local.upper[d] + gw[d];
			}
			else if(gc && d < di)
			{
				ghost_back [di].lower[d] = local.lower[d] - gw[d];
				ghost_back [di].upper[d] = local.upper[d] + gw[d];
				ghost_front[di].lower[d] = local.lower[d] - gw[d];
				ghost_front[di].upper[d] = local.upper[d] + gw[d];
			}
			else
			{
				ghost_back [di].lower[d] = local.lower[d];
				ghost_back [di].upper[d] = local.upper[d];
				ghost_front[di].lower[d] = local.lower[d];
				ghost_front[di].upper[d] = local.upper[d];
			}

#ifndef NDEBUG
	if(log_mask[verbose_alloc])
		log_out() << "allocate " << local << endl;
#endif
#if defined(SHARK_MPI_COMM)
	{
		// Create ghost types
		MPI_Datatype base = mpi_type<T>::t;
		if(mpi_type<T>::count() != 1)
			MPI_Type_contiguous(mpi_type<T>::count(), base, &base);
		MPI_Aint lb, extent;
		MPI_Type_get_extent(base, &lb, &extent);

		for(int di = 0; di < ndim; di++)
		{
			MPI_Type_dup(base, &impl->ghost[di]);

			for(int d = ndim-1; d >= 0; d--)
			{
				MPI_Datatype tmp = impl->ghost[di];
				coord n = ghost_back[di].upper[d] - ghost_back[di].lower[d];
				MPI_Type_create_hvector(n, 1, ld[d+1]*extent, tmp, &impl->ghost[di]);
				MPI_Type_free(&tmp);
			}

			MPI_Type_commit(&impl->ghost[di]);
		}

		if(mpi_type<T>::count() != 1)
			MPI_Type_free(&base);
	}

	// Allocate memory
	MPI_Aint size = ld[0] * sizeof(T);
	MPI_Alloc_mem(size, MPI_INFO_NULL, &ptr);

	// Create window
	MPI_Win_create(ptr, size, static_cast<int>(sizeof(T)), MPI_INFO_NULL, domain().group.impl->comm, &impl->win);
#elif defined(SHARK_NO_COMM)
	unused(impl);
	ptr = static_cast<T*>(mem_alloc(ld[0] * sizeof(T)));
#else
#error "No comm allocate"
#endif

	{
		typename Domain<ndim>::pcoords np = domain().np;
		typename Domain<ndim>::pcoords ip = domain().indexp();
		Access<ndim,T> acc(*this);

		for(int di = 0; di < ndim; di++)
		{
			typename Boundary<ndim,T>::fixed_type* b = dynamic_cast<typename Boundary<ndim,T>::fixed_type*>(bd[di].t.get());

			if(b != nullptr)
			{
				if(ip[di] == 0)
					b->set(acc, ghost_back[di]);
				if(ip[di] == np[di]-1)
					b->set(acc, ghost_front[di]);
			}
		}
	}
}

// Constructors
template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(const Domain<ndim>& domain, coords<ndim> ghost_width, bool ghost_corners, bounds bd):
  dom(&domain), gw(ghost_width), gc(ghost_corners), bd(bd), impl(new GAImpl<ndim>()), lc(0) {
	allocate();
}

template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(const GlobalArray<ndim,T>& other, bool copy):
  dom(other.dom), gw(other.gw), gc(other.gc), bd(other.bd), impl(new GAImpl<ndim>()), lc(0) {
	allocate();

	if(copy)
		*this = other;
}

// Destructor

template<int ndim,typename T>
void GlobalArray<ndim,T>::deallocate()
{
#if defined(SHARK_MPI_COMM)
	for(int di = 0; di < ndim; di++)
		MPI_Type_free(&impl->ghost[di]);

	MPI_Win_free(&impl->win);
	MPI_Free_mem(ptr);
#elif defined(SHARK_NO_COMM)
	mem_free(ptr);
#else
#error "No comm deallocate"
#endif
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::reshape(const Domain<ndim>& domain)
{
	assert(this->domain().equiv(domain));

	GlobalArray<ndim,T> tmp(domain, gw, gc, bd);

	domain.group.sync();
	if(domain.hasData())
		this->get(domain.local(), essential_lead<ndim>(tmp.ld), tmp.ptr + gw.offset(tmp.ld));
	domain.group.sync();

	*this = move(tmp);
}

#ifdef SHARK_MPI_COMM

namespace {
	template<int ndim>
	class UpdateHandle: public Handle {
		template<int,typename> friend class GlobalArray;
		MPI_Request req[4*ndim];
	public:
		UpdateHandle();
		virtual ~UpdateHandle();
		virtual bool test();
		virtual void wait();
	};

	template<int ndim, typename T>
	class GatherHandle: public Handle {
		friend class GlobalArray<ndim,T>;
		vector<MPI_Request> reqs;
		const Access<ndim,T> acc;
		bool testReqs();
	public:
		GatherHandle(const GlobalArray<ndim,T>& ga);
		virtual ~GatherHandle();
		virtual bool test();
		virtual void wait();
	};


	template<int ndim, typename T>
	class ScatterHandle: public Handle
	{
		friend class GlobalArray<ndim,T>;
		GlobalArray<ndim,T>& ga;
		vector<MPI_Request> reqs_local;
		vector<MPI_Request> reqs_global;
		vector<pair<coords_range<ndim>,unique_ptr<T[]>>> bufs;
		void accumulate(int i);
		bool testLocal();
		bool testGlobal();
	public:
		ScatterHandle(GlobalArray<ndim,T>& ga);
		virtual ~ScatterHandle();
		virtual bool test();
		virtual void wait();
	};
}

template<int ndim>
UpdateHandle<ndim>::UpdateHandle()
{
	for(int di = 0; di < ndim; di++)
	{
		req[4*di] = MPI_REQUEST_NULL;
		req[4*di+1] = MPI_REQUEST_NULL;
		req[4*di+2] = MPI_REQUEST_NULL;
		req[4*di+3] = MPI_REQUEST_NULL;
	}
}

template<int ndim>
UpdateHandle<ndim>::~UpdateHandle()
{
}

template<int ndim>
void UpdateHandle<ndim>::wait()
{
	MPI_Waitall(4*ndim, req, MPI_STATUSES_IGNORE);
}

template<int ndim>
bool UpdateHandle<ndim>::test()
{
	int flag;
	MPI_Testall(4*ndim, req, &flag, MPI_STATUSES_IGNORE);
	return flag;
}

template<int ndim,typename T>
GatherHandle<ndim,T>::GatherHandle(const GlobalArray<ndim,T>& ga): acc(ga) {
}

template<int ndim,typename T>
GatherHandle<ndim,T>::~GatherHandle() {
}

template<int ndim,typename T>
void GatherHandle<ndim,T>::wait() {
	if(!reqs.empty())
		MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
}

template<int ndim,typename T>
bool GatherHandle<ndim,T>::testReqs() {
	int flag;
	MPI_Testall(reqs.size(), reqs.data(), &flag, MPI_STATUSES_IGNORE);
	return flag;
}

template<int ndim,typename T>
bool GatherHandle<ndim,T>::test() {
	return reqs.empty() || testReqs();
}

template<int ndim,typename T>
ScatterHandle<ndim,T>::ScatterHandle(GlobalArray<ndim,T>& ga): ga(ga) {
}

template<int ndim,typename T>
ScatterHandle<ndim,T>::~ScatterHandle() {
}

template<int ndim,typename T>
void ScatterHandle<ndim,T>::accumulate(int i) {
	Access<ndim,T> acc(ga);
	const coords_range<ndim>& r(bufs[i].first);
	const T* buf(bufs[i].second.get());
	const coords<ndim+1> ld(r.stride());
	r.for_each([&acc,&r,&buf,&ld](coords<ndim> ii) {
		acc(ii) += buf[(ii - r.lower).offset(ld)];
	});
}

template<int ndim,typename T>
void ScatterHandle<ndim,T>::wait() {
	if(!reqs_global.empty())
		while(true) {
			int i;
			MPI_Waitany(reqs_global.size(), reqs_global.data(), &i, MPI_STATUS_IGNORE);
			if(i == MPI_UNDEFINED)
				break;
			accumulate(i);
		}
	if(!reqs_local.empty())
		MPI_Waitall(reqs_local.size(), reqs_local.data(), MPI_STATUSES_IGNORE);
}

template<int ndim,typename T>
bool ScatterHandle<ndim,T>::testLocal() {
	int flag;
	MPI_Testall(reqs_local.size(), reqs_local.data(), &flag, MPI_STATUSES_IGNORE);
	return flag;
}

template<int ndim,typename T>
bool ScatterHandle<ndim,T>::testGlobal() {
	int flag;
	// Try to make progress
	while(true) {
		int i;
		MPI_Testany(reqs_global.size(), reqs_global.data(), &i, &flag, MPI_STATUS_IGNORE);
		if(i == MPI_UNDEFINED)
			break;
		accumulate(i);
	}
	// No more progress
	// flag = true: all requests inactive
	// flag = false: still active, incomplete requests left
	return flag;
}

template<int ndim,typename T>
bool ScatterHandle<ndim,T>::test() {
	return (reqs_global.empty() || testGlobal()) && (reqs_local.empty() || testLocal());
}

#endif

template<int ndim,typename T>
void GlobalArray<ndim,T>::update(long k) const {
	iupdate(k).wait();
}

template<int ndim,typename T>
Future<void> GlobalArray<ndim,T>::iupdate(long k) const
{
	const coords<ndim> gw = ghost_width();

	typename Domain<ndim>::pcoords np = domain().np;
	typename Domain<ndim>::pcoords ip = domain().indexp();

	Access<ndim,T> acc(*this);
	
#if defined(SHARK_MPI_COMM)
	const bool gc = ghost_corners();
	MPI_Comm comm = domain().group.impl->comm;
	unique_ptr<UpdateHandle<ndim>> h(new UpdateHandle<ndim>());
#endif

	for(int di = 0; di < ndim; di++)
	{
		if(gw[di] > 0)
		{
#ifndef NDEBUG
			if(log_mask[verbose_update])
			{
				log_out() << "update back: " << ghost_back[di] << endl;
				log_out() << "update front: " << ghost_front[di] << endl;
			}
#endif
			// 1) Consider updating ghosts based on boundary
			typename Boundary<ndim,T>::general_type* b = dynamic_cast<typename Boundary<ndim,T>::general_type*>(bd[di].t.get());

			if(b != nullptr)
			{
				if(ip[di] == 0)
					b->set(acc, ghost_back[di], k);
				if(ip[di] == np[di]-1)
					b->set(acc, ghost_front[di], k);
			}

			// 2) Consider updating ghosts based on communication with neighbors (or periodic equivalents)
			bool pd = dynamic_cast<typename Boundary<ndim,T>::periodic_type*>(bd[di].t.get()) != nullptr;
#if defined(SHARK_MPI_COMM)
			int prev = domain().shiftd(di, -1, pd);
			int next = domain().shiftd(di,  1, pd);
			// backward
			MPI_Isend(&acc(ghost_back[di].adj(di, 1)), 1, impl->ghost[di], prev, 2*di, comm, &h->req[4*di]);
			MPI_Irecv(&acc(ghost_front[di].lower), 1, impl->ghost[di], next, 2*di, comm, &h->req[4*di+1]);
			// forward
			MPI_Isend(&acc(ghost_front[di].adj(di, -1)), 1, impl->ghost[di], next, 2*di+1, comm, &h->req[4*di+2]);
			MPI_Irecv(&acc(ghost_back[di].lower), 1, impl->ghost[di], prev, 2*di+1, comm, &h->req[4*di+3]);
			if(gc)
			{
				std::cerr << " WAITING " <<std::endl;
				MPI_Waitall(4, &h->req[4*di], MPI_STATUSES_IGNORE);
			}
#elif defined(SHARK_NO_COMM)
			if(pd) {
				{
					coords<ndim> off = ghost_back[di].adj(di, 1) - ghost_front[di].lower;
					ghost_front[di].for_each([&acc,&off](coords<ndim> ii) {
						acc(ii) = acc(ii+off);
					});
				}
				{
					coords<ndim> off = ghost_front[di].adj(di, -1) - ghost_back[di].lower;
					ghost_back[di].for_each([&acc,&off](coords<ndim> ii) {
						acc(ii) = acc(ii+off);
					});
				}
			}
#else
#error "No comm update"
#endif

		}
	}

#if defined(SHARK_MPI_COMM)
	if(!gc)
		return Future<void>(std::move(h));
	else
		return Future<void>(make_unique<DoneHandle>());
#elif defined(SHARK_NO_COMM)
	return Future<void>(make_unique<DoneHandle>());
#endif
}

//IMEN BEGIN

template<int ndim,typename T>
MPI_Request* GlobalArray<ndim,T>::updateBegin(long k) const
{
	const coords<ndim> gw = ghost_width();

	typename Domain<ndim>::pcoords np = domain().np;
	typename Domain<ndim>::pcoords ip = domain().indexp();

	Access<ndim,T> acc(*this);

	#if defined(SHARK_MPI_COMM)
		MPI_Comm comm = domain().group.impl->comm;
		unique_ptr<UpdateHandle<ndim>> h(new UpdateHandle<ndim>());
	#endif

		for(int di = 0; di < ndim; di++)
		{
			if(gw[di] > 0)
			{
	#ifndef NDEBUG
				if(log_mask[verbose_update])
				{
					log_out() << "update back: " << ghost_back[di] << endl;
					log_out() << "update front: " << ghost_front[di] << endl;
				}
	#endif
				// 1) Consider updating ghosts based on boundary
				typename Boundary<ndim,T>::general_type* b = dynamic_cast<typename Boundary<ndim,T>::general_type*>(bd[di].t.get());

				if(b != nullptr)
				{
					if(ip[di] == 0)
						b->set(acc, ghost_back[di], k);
					if(ip[di] == np[di]-1)
						b->set(acc, ghost_front[di], k);
				}

				// 2) Consider updating ghosts based on communication with neighbors (or periodic equivalents)
				bool pd = dynamic_cast<typename Boundary<ndim,T>::periodic_type*>(bd[di].t.get()) != nullptr;
	#if defined(SHARK_MPI_COMM)
				int prev = domain().shiftd(di, -1, pd);
				int next = domain().shiftd(di,  1, pd);
				// backward
				MPI_Isend(&acc(ghost_back[di].adj(di, 1)), 1, impl->ghost[di], prev, 2*di, comm, &h->req[4*di]);
				MPI_Irecv(&acc(ghost_front[di].lower), 1, impl->ghost[di], next, 2*di, comm, &h->req[4*di+1]);
				// forward
				MPI_Isend(&acc(ghost_front[di].adj(di, -1)), 1, impl->ghost[di], next, 2*di+1, comm, &h->req[4*di+2]);
				MPI_Irecv(&acc(ghost_back[di].lower), 1, impl->ghost[di], prev, 2*di+1, comm, &h->req[4*di+3]);

	#elif defined(SHARK_NO_COMM)
				if(pd)
				{
					{
						coords<ndim> off = ghost_back[di].adj(di, 1) - ghost_front[di].lower;
						ghost_front[di].for_each([&acc,&off](coords<ndim> ii) {
							acc(ii) = acc(ii+off);
						});
					}
					{
						coords<ndim> off = ghost_front[di].adj(di, -1) - ghost_back[di].lower;
						ghost_back[di].for_each([&acc,&off](coords<ndim> ii) {
							acc(ii) = acc(ii+off);
						});
					}
				}
	#else
	#error "No comm update"
	#endif

			}
		}

		return h->req;
}

template<int ndim,typename T>
void GlobalArray<ndim,T>::updateWait(MPI_Request* req) const
{
	MPI_Waitall(4*ndim, req, MPI_STATUSES_IGNORE);
}

//IMEN END

template<int ndim, typename T>
void GlobalArray<ndim,T>::get(coords_range<ndim> range, T* buf) const {
	get(range, essential_lead<ndim>(range.stride()), buf);
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::get(coords_range<ndim> range, array<size_t,ndim-1> ld, T* buf) const {
#if defined(SHARK_MPI_COMM)
	typename Domain<ndim>::ProcessOverlap(domain(), range).visit([this,range,ld,buf](int id, coords_range<ndim> i) {
#ifndef NDEBUG
			if(log_mask[verbose_rma])
				this->log_out() << "get " << id << ": " << i << endl;
#endif
			// Target
			coords_range<ndim> tgt = this->domain().local(id);
			coords<ndim+1> tgtld = tgt.stride(gw);
			coord tgtoff = (i.lower - tgt.lower + gw).offset(tgtld);
			mpi_type_block<ndim,T> tgtt(i.counts(), tgtld);
			// Origin
			coord ognoff = (i.lower - range.lower).offset(ld);
			mpi_type_block<ndim,T> ognt(i.counts(), ld);
			// RMA
			MPI_Win_lock(MPI_LOCK_SHARED, id, 0, impl->win);
			MPI_Get(buf+ognoff, 1, ognt.t, id, tgtoff, 1, tgtt.t, impl->win);
			MPI_Win_unlock(id, impl->win);
	});
#elif defined(SHARK_NO_COMM)
#ifndef NDEBUG
	if(log_mask[verbose_rma])
		log_out() << "get " << range << endl;
#endif
	const Access<ndim,T> acc(*this);
	range.for_each([&acc,&buf,&range,&ld](coords<ndim> i) {
		buf[(i - range.lower).offset(ld)] = acc(i);
	});
#else
#error "No comm get"
#endif
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::put(coords_range<ndim> range, const T* buf) {
	put(range, essential_lead<ndim>(range.stride()), buf);
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::put(coords_range<ndim> range, array<size_t,ndim-1> ld, const T* buf) {
#if defined(SHARK_MPI_COMM)
	typename Domain<ndim>::ProcessOverlap(domain(), range).visit([this,range,ld,buf](int id, coords_range<ndim> i) {
#ifndef NDEBUG
			if(log_mask[verbose_rma])
				this->log_out() << "put " << id << ": " << i << endl;
#endif
			// Target
			coords_range<ndim> tgt = this->domain().local(id);
			coords<ndim+1> tgtld = tgt.stride(gw);
			coord tgtoff = (i.lower - tgt.lower + gw).offset(tgtld);
			mpi_type_block<ndim,T> tgtt(i.counts(), tgtld);
			// Origin
			coord ognoff = (i.lower - range.lower).offset(ld);
			mpi_type_block<ndim,T> ognt(i.counts(), ld);
			// RMA
			MPI_Win_lock(MPI_LOCK_EXCLUSIVE, id, 0, impl->win);
			MPI_Put(const_cast<T*>(buf+ognoff), 1, ognt.t, id, tgtoff, 1, tgtt.t, impl->win);
			MPI_Win_unlock(id, impl->win);
	});
#elif defined(SHARK_NO_COMM)
#ifndef NDEBUG
	if(log_mask[verbose_rma])
		log_out() << "put " << range << endl;
#endif
	Access<ndim,T> acc(*this);
	range.for_each([&acc,&buf,&range,&ld](coords<ndim> i) {
		acc(i) = buf[(i - range.lower).offset(ld)];
	});
#else
#error "No comm put"
#endif
}

template<int ndim, typename T> template<typename>
void GlobalArray<ndim,T>::accumulate(coords_range<ndim> range, const T* buf) {
	accumulate(range, essential_lead<ndim>(range.stride()), buf);
}

template<int ndim,typename T> template<typename>
void GlobalArray<ndim,T>::accumulate(coords_range<ndim> range, array<size_t,ndim-1> ld, const T* buf) {
#if defined(SHARK_MPI_COMM)
	typename Domain<ndim>::ProcessOverlap(domain(), range).visit([this,range,ld,buf](int id, coords_range<ndim> i) {

#if defined(OMPI_MAJOR_VERSION) && (OMPI_MAJOR_VERSION < 2) && (OMPI_MINOR_VERSION < 8)
#error MPI_Accumulate does not seem to work with vector types on OpenMPI < 1.8
#endif

#ifndef NDEBUG
			if(log_mask[verbose_rma])
				this->log_out() << "accumulate " << id << ": " << i << endl;
#endif
			// Target
			coords_range<ndim> tgt = this->domain().local(id);
			coords<ndim+1> tgtld = tgt.stride(gw);
			coord tgtoff = (i.lower - tgt.lower + gw).offset(tgtld);
			mpi_type_block<ndim,T> tgtt(i.counts(), tgtld);
			// Origin
			coord ognoff = (i.lower - range.lower).offset(ld);
			mpi_type_block<ndim,T> ognt(i.counts(), ld);
			// RMA
			MPI_Win_lock(MPI_LOCK_SHARED, id, 0, impl->win);
			MPI_Accumulate(const_cast<T*>(buf+ognoff), 1, ognt.t, id, tgtoff, 1, tgtt.t, MPI_SUM, impl->win);
			MPI_Win_unlock(id, impl->win);
	});
#elif defined(SHARK_NO_COMM)
#ifndef NDEBUG
	if(log_mask[verbose_rma])
		log_out() << "accumulate " << range << endl;
#endif
	Access<ndim,T> acc(*this);
	range.for_each([&acc,&buf,&range,&ld](coords<ndim> i) {
		acc(i) += buf[(i - range.lower).offset(ld)];
	});
#else
#error "No comm accumulate"
#endif
}

template<int ndim,typename T> template<typename>
Future<void> GlobalArray<ndim,T>::igather(SparseArray<ndim,T>& sa) const {
	assert(domain() == sa.dom);
#if defined(SHARK_NO_COMM)
	auto eld = essential_lead<ndim>(sa.ld);
	domain().sync();
	sa.iter([this,&sa,&eld](const coords_range<ndim>& r) {
		this->get(r, eld, sa.ptr + r.lower.offset(sa.ld));
	});
	domain().sync();
	return Future<void>(make_unique<DoneHandle>());
#elif defined(SHARK_MPI_COMM)
	const Domain<ndim>& dom(domain());
	int nprocs = dom.group.impl->size();
	MPI_Comm comm = dom.group.impl->comm;
	vector<vector<coords_range<ndim>>> local_ranges(nprocs), global_ranges(nprocs);
	sa.proc_ranges(local_ranges, global_ranges);
#ifndef NDEBUG
	if(log_mask[verbose_collective])
		for(int k = 0; k < nprocs; k++)
			for(auto it = local_ranges[k].cbegin(); it != local_ranges[k].cend(); ++it)
				log_out() << "gather from " << k << ": " << *it << endl;
#endif
	unique_ptr<GatherHandle<ndim,T>> h(new GatherHandle<ndim,T>(*this));
	for(int k = 0; k < nprocs; k++)
		for(auto it = local_ranges[k].cbegin(); it != local_ranges[k].cend(); ++it) {
			h->reqs.emplace_back();
			MPI_Irecv(&sa.ptr[it->lower.offset(sa.ld)], 1, mpi_type_block<ndim,T>(it->counts(), sa.ld).t, k, 0, comm, &h->reqs.back());
		}
	for(int k = 0; k < nprocs; k++)
		for(auto it = global_ranges[k].cbegin(); it != global_ranges[k].cend(); ++it) {
			h->reqs.emplace_back();
			MPI_Isend(const_cast<T*>(&h->acc(it->lower)), 1, mpi_type_block<ndim,T>(it->counts(), ld).t, k, 0, comm, &h->reqs.back());
		}
	return Future<void>(std::move(h));
#else
#error "No comm gather"
#endif
}

template<int ndim,typename T> template<typename>
void GlobalArray<ndim,T>::gather(SparseArray<ndim,T>& sa) const {
	igather(sa).wait();
}

template<int ndim,typename T> template<typename>
Future<void> GlobalArray<ndim,T>::iscatterAcc(const SparseArray<ndim,T>& sa) {
	assert(domain() == sa.dom);
#if defined(SHARK_NO_COMM)
	auto eld = essential_lead<ndim>(sa.ld);
	domain().sync();
	sa.iter([this,&sa,&eld](const coords_range<ndim>& r) {
		this->accumulate(r, eld, sa.ptr + r.lower.offset(sa.ld));
	});
	domain().sync();
	return Future<void>(make_unique<DoneHandle>());
#elif defined(SHARK_MPI_COMM)
	const Domain<ndim>& dom(domain());
	int nprocs = dom.group.impl->size();
	MPI_Comm comm = dom.group.impl->comm;
	vector<vector<coords_range<ndim>>> local_ranges(nprocs), global_ranges(nprocs);
	sa.proc_ranges(local_ranges, global_ranges);
#ifndef NDEBUG
	if(log_mask[verbose_collective])
		for(int k = 0; k < nprocs; k++)
			for(auto it = local_ranges[k].cbegin(); it != local_ranges[k].cend(); ++it)
				log_out() << "scatterAcc to " << k << ": " << *it << endl;
#endif
	unique_ptr<ScatterHandle<ndim,T>> h(new ScatterHandle<ndim,T>(*this));
	for(int k = 0; k < nprocs; k++)
		for(auto it = global_ranges[k].cbegin(); it != global_ranges[k].cend(); ++it) {
			h->bufs.emplace_back(*it, unique_ptr<T[]>(new T[it->count()]));
			h->reqs_global.emplace_back();
			MPI_Irecv(h->bufs.back().second.get(), it->count() * mpi_type<T>::count(), mpi_type<T>::t, k, 0, comm, &h->reqs_global.back());
		}
	for(int k = 0; k < nprocs; k++)
		for(auto it = local_ranges[k].cbegin(); it != local_ranges[k].cend(); ++it) {
			h->reqs_local.emplace_back();
			MPI_Isend(&sa.ptr[it->lower.offset(sa.ld)], 1, mpi_type_block<ndim,T>(it->counts(), sa.ld).t, k, 0, comm, &h->reqs_local.back());
		}
	return Future<void>(std::move(h));
#else
#error "No comm scatterAcc"
#endif
}
template<int ndim,typename T> template<typename>
void GlobalArray<ndim,T>::dump(coords_range<ndim> range, std::string filename) {
#if defined(SHARK_MPI_COMM)
#ifndef NDEBUG
        if(log_mask[verbose_rma])
                this->log_out() << "dump to " << filename << endl;
#endif
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

        // type for part to read from
        coords_range<ndim> src = domain().local().overlap(range);
        coords<ndim+1> srcld = src.stride(gw);
        const T* buf = &da(src.lower);
        mpi_type_block<ndim,T> srct(src.counts(), srcld);
        
        // where to write to
        coords_range<ndim> dst = region().overlap(range);
        coords<ndim+1> dstld = dst.stride(); // no ghost regions in file
        coord dstoff = dst.lower.offset(dstld); // offset wrt to global array
        
        MPI_Status stat;
        MPI_File_write_at(fh, dstoff, buf, 1, srct.t, &stat);

        MPI_File_close(&fh);
#elif defined(SHARK_NO_COMM)
#warning "Dump not implemented for SHARK_NO_COMM
#else
#error "No comm dump"
#endif
}

template<int ndim,typename T> template<typename>
void GlobalArray<ndim,T>::scatterAcc(const SparseArray<ndim,T>& sa) {
	iscatterAcc(sa).wait();
}

template<int ndim,typename T>
GADest<ndim,T> GlobalArray<ndim,T>::region(coords_range<ndim> r) {
	return GADest<ndim,T>(*this, r);
}

template<int ndim,typename T>
GADest<ndim,T>::GADest(GlobalArray<ndim,T>& ga, coords_range<ndim> r): ga(ga), r(r) { }

template<int ndim,typename T>
GADest<ndim,T>::~GADest() { }

template<int ndim,typename T>
GADest<ndim,T>::GADest(const GADest& gad): ga(gad.ga), r(gad.r) { }

template<int ndim, typename T>
GARef<ndim,T>::GARef(const GlobalArray<ndim,T>& ga): ga(ga) { }

template<int ndim, typename T>
GARef<ndim,T>::~GARef() { }

// Set-up instantiations

#include "types"

#define SYMBDT(d,T) template class shark::ndim::GlobalArray<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

#define SYMB_ARITH

#define SYMBDT(d,T) template void shark::ndim::GlobalArray<d,T>::accumulate(coords_range<d>, const T*);
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template void shark::ndim::GlobalArray<d,T>::accumulate(coords_range<d>, std::array<std::size_t,d-1>, const T*);
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template void shark::ndim::GlobalArray<d,T>::gather(SparseArray<d,T>&) const;
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template Future<void> shark::ndim::GlobalArray<d,T>::igather(SparseArray<d,T>&) const;
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template void shark::ndim::GlobalArray<d,T>::scatterAcc(const SparseArray<d,T>&);
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template Future<void> shark::ndim::GlobalArray<d,T>::iscatterAcc(const SparseArray<d,T>&);
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template void  shark::ndim::GlobalArray<d,T>::dump(coords_range<d> range, std::string);
#include "inst_dimtype"
#undef SYMBDT

#undef SYMB_ARITH

#define SYMBDT(d,T) template class shark::ndim::GADest<d,T>; 
#include "inst_dimtype"
#undef SYMBDT
	
#define SYMBDT(d,T) template class shark::ndim::GARef<d,T>; 
#include "inst_dimtype"
#undef SYMBDT
