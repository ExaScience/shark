/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <set>
#include <cstring>
#include <map>

#include <shark/access.hpp>
#include <shark/boundary.hpp>
#include <shark/domain.hpp>
#include <shark/globalarray.hpp>
#include <shark/globals.hpp>
#include <shark/common.hpp>

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
	return *shark::log_out << "GA:[" << ptr << "]"
#if defined(SHARK_GPI_COMM)
            << (int)impl->seg << "@"
#endif
            << domain().group.procid << ":";


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
	lc(0),
	g_out_back(std::move(other.g_out_back)),
	g_out_front(std::move(other.g_out_front)),
	g_in_back(std::move(other.g_in_back)),
	g_in_front(std::move(other.g_in_front))
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
	g_out_back = std::move(other.g_out_back);
	g_out_front = std::move(other.g_out_front);
	g_in_back = std::move(other.g_in_back);
	g_in_front = std::move(other.g_in_front);
	other.reset();
	return *this;
}

template<int ndim, typename T>
GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator=(const GlobalArray<ndim,T>& other) {
	region(domain().total(), false) = other;
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

	// Allocate memory
	ld = local.stride(gw);
        unsigned count = ld[0];

	for(int di = 0; di < ndim; di++) {
		for(int d = 0; d < ndim; d++) {
#ifdef SHARK_STPACK
                        coord gw_prev = gw[d];
                        coord gw_next = gw[d];
#else
                        bool pd = dynamic_cast<typename Boundary<ndim,T>::periodic_type*>(bd[di].t.get()) != nullptr;
                        int prev = domain().shiftd(di, -1, pd);
                        int next = domain().shiftd(di,  1, pd);
                        coord gw_prev = prev < 0 ? 0 : gw[d];
                        coord gw_next = next < 0 ? 0 : gw[d];
#endif

			if(d == di)
			{
				ghost_back [di].lower[d] = local.lower[d] - gw_prev;
				ghost_back [di].upper[d] = local.lower[d];
				ghost_front[di].lower[d] = local.upper[d];
				ghost_front[di].upper[d] = local.upper[d] + gw_next;
			}
			else if(gc && d < di)
			{
				ghost_back [di].lower[d] = local.lower[d] - gw_prev;
				ghost_back [di].upper[d] = local.upper[d] + gw_prev;
				ghost_front[di].lower[d] = local.lower[d] - gw_next;
				ghost_front[di].upper[d] = local.upper[d] + gw_next;
			}
			else
			{
				ghost_back [di].lower[d] = local.lower[d];
				ghost_back [di].upper[d] = local.upper[d];
				ghost_front[di].lower[d] = local.lower[d];
				ghost_front[di].upper[d] = local.upper[d];
			}
                }

                count += 2* ghost_back[di].count();
                count += 2* ghost_front[di].count();
        }

	unsigned size = count * sizeof(T);

#ifndef NDEBUG
	if(log_mask[verbose_alloc])
		log_out() << "allocate " << local << endl;
#endif
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
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


#if defined(SHARK_GPI_COMM)
        impl->seg = GroupImpl::next_seg_id();
        success_or_die(gaspi_segment_create(impl->seg, size, domain().group.impl->grp, GASPI_BLOCK, GASPI_ALLOC_DEFAULT));  
        gaspi_segment_ptr(impl->seg, (gaspi_pointer_t*)&ptr);
        assert(ptr);
#else
	MPI_Alloc_mem(size, MPI_INFO_NULL, &ptr);
#endif

	// Create window
	MPI_Win_create(ptr, size, static_cast<int>(sizeof(T)), MPI_INFO_NULL, domain().group.impl->comm, &impl->win);
#elif defined(SHARK_NO_COMM)
	unused(impl);
	ptr = static_cast<T*>(mem_alloc(size));
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

                    if(b == nullptr) continue;
                    if(ip[di] == 0) b->set(acc, ghost_back[di]);
                    if(ip[di] == np[di]-1) b->set(acc, ghost_front[di]);
                }
	}

        auto p = ptr + ld[0];
	for(int di = 0; di < ndim; di++) {
                g_in_back  [di] = GABuf<ndim,T>(ghost_in_back(di), p); 
                p += ghost_in_back(di).count();
                g_in_front [di] = GABuf<ndim,T>(ghost_in_front(di), p);
                p += ghost_in_front(di).count();
                g_out_back [di] = GABuf<ndim,T>(ghost_out_back(di), p);
                p += ghost_out_back(di).count();
                g_out_front[di] = GABuf<ndim,T>(ghost_out_front(di), p);
                p += ghost_out_front(di).count();
        }
        assert(p == ptr + count);
}

// Constructors
template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(const Domain<ndim>& domain, coords<ndim> ghost_width, bool ghost_corners, bounds bd):
  dom(&domain), gw(ghost_width), gc(ghost_corners), bd(bd), impl(new GAImpl<ndim,T>()), lc(0) {
	allocate();
}

template<int ndim, typename T>
GlobalArray<ndim,T>::GlobalArray(const GlobalArray<ndim,T>& other, bool copy):
  dom(other.dom), gw(other.gw), gc(other.gc), bd(other.bd), impl(new GAImpl<ndim,T>()), lc(0) {
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
#elif defined(SHARK_GPI_COMM)
	gaspi_segment_delete(impl->seg);
        GroupImpl::free_seg_id(impl->seg);
#elif defined(SHARK_NO_COMM)
	mem_free(ptr);
#else
#error "No comm deallocate"
#endif
}

#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)

namespace {
	template<int ndim>
	class UpdateHandle: public Handle {
	public:

		template<int,typename> friend class GlobalArray;
		MPI_Request req[4*ndim];

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
        SHARK_COUNTER("test");
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

#if defined(SHARK_GPI_COMM)

namespace {
	template<int ndim>
	class GPIHandle: public Handle {
        private:
                const gaspi_segment_id_t seg;
                const gaspi_notification_t expected;
                std::set<gaspi_notification_id_t> ids;

                int from() const { return *ids.cbegin(); }
                int to()   const { return *ids.crbegin(); } 
                int range() const { return to() - from() + 1; }

	public:
		GPIHandle(gaspi_segment_id_t, gaspi_notification_t);
		virtual ~GPIHandle();
		virtual bool test();
		virtual void wait();
                GPIHandle *add(gaspi_notification_id_t id) { ids.insert(id); return this; }
	};
}

template<int ndim>
GPIHandle<ndim>::GPIHandle(gaspi_segment_id_t s, gaspi_notification_t e)
    : seg(s), expected(e) { }

template<int ndim>
GPIHandle<ndim>::~GPIHandle()
{
        assert(ids.empty());

        gaspi_notification_t v;
        for ( auto i : ids ) {
                SUCCESS_OR_DIE(gaspi_notify_reset(seg, i, &v));
        }
}

template<int ndim>
void GPIHandle<ndim>::wait()
{
        while (!ids.empty()) {
            gaspi_notification_id_t id;
            SUCCESS_OR_DIE(gaspi_notify_waitsome(seg, from(), range(), &id, GASPI_BLOCK));

            if (ids.find(id) == ids.end()) continue;
            ids.erase(id);
            gaspi_notification_t v;
            SUCCESS_OR_DIE(gaspi_notify_reset(seg, id, &v));
            assert(v == expected);
        }
}

template<int ndim>
bool GPIHandle<ndim>::test()
{
        for ( auto i : ids ) {
            gaspi_notification_id_t id;
            auto r = gaspi_notify_waitsome(seg, i, 1, &id, 0);
            if (r != GASPI_SUCCESS) return false;
            ids.erase(id);
            gaspi_notification_t v;
            SUCCESS_OR_DIE(gaspi_notify_reset(seg, id, &v));
            assert(v == expected);
        }

        return true;
}

#endif

///---------------------------------------------------

template<int ndim,typename T>
void GlobalArray<ndim,T>::update(long k) const {
    SHARK_COUNTER("update");
	iupdate(k).wait();
}


template<int ndim,typename T>
std::ostream& shark::ndim::operator<<(std::ostream& os, const GlobalArray<ndim,T> &ga) {
    Access<ndim,T> acc(ga);
    os << ga.domain().local() << ": [";
    ga.domain().local().for_each([&os,&acc](coords<ndim> ii) { os << acc(ii) << ", "; });
    os << "]";
    return os;
}

template<int ndim,typename T>
Future<void> GlobalArray<ndim,T>::iupdate(long k) const
{
        SHARK_COUNTER("iupdate");
	const coords<ndim> gw = ghost_width();

	typename Domain<ndim>::pcoords np = domain().np;
	typename Domain<ndim>::pcoords ip = domain().indexp();

	Access<ndim,T> acc(*this);
	
#if defined(SHARK_MPI_COMM)
	const bool gc = ghost_corners();
	MPI_Comm comm = domain().group.impl->comm;
	unique_ptr<UpdateHandle<ndim>> h(new UpdateHandle<ndim>());
#endif

#if defined(SHARK_GPI_COMM)
        static int notify_value = 0;
        auto handle = new GPIHandle<ndim>(impl->seg, ++notify_value);
	unique_ptr<GPIHandle<ndim>> h(handle);
#endif

	for(int di = 0; di < ndim; di++) {
		if(gw[di] <= 0) continue;
#ifndef NDEBUG
                if(log_mask[verbose_update]) {
                        Access<ndim,T> acc(*this);
                        log_out() << "update back: " << ghost_back[di] << endl;
                        log_out() << "update front: " << ghost_front[di] << endl;
                }
#endif
                // 1) Consider updating ghosts based on boundary
                typename Boundary<ndim,T>::general_type* b = dynamic_cast<typename Boundary<ndim,T>::general_type*>(bd[di].t.get());
                if(b != nullptr) {
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
#ifdef SHARK_STPACK
                // backward
                MPI_Isend(&acc(ghost_out_back(di).lower), 1, impl->ghost[di], prev, 2*di, comm, &h->req[4*di]);
                MPI_Irecv(&acc(ghost_in_front(di).lower), 1, impl->ghost[di], next, 2*di, comm, &h->req[4*di+1]);
                // forward
                MPI_Isend(&acc(ghost_out_front(di).lower), 1, impl->ghost[di], next, 2*di+1, comm, &h->req[4*di+2]);
                MPI_Irecv(&acc(ghost_in_back(di).lower), 1, impl->ghost[di], prev, 2*di+1, comm, &h->req[4*di+3]);
#else
                // backward
                const int size_out_back   = g_out_back[di].r.count();
                const int size_out_front  = g_out_front[di].r.count();

                MPI_Isend(g_out_back[di].ptr, size_out_back, mpi_type<T>::t, prev, 2*di, comm, &h->req[4*di]);
                MPI_Irecv(g_in_front[di].ptr, size_out_front, mpi_type<T>::t, next, 2*di, comm, &h->req[4*di+1]);
                // forward
                MPI_Isend(g_out_front[di].ptr, size_out_front, mpi_type<T>::t, next, 2*di+1, comm, &h->req[4*di+2]);
                MPI_Irecv(g_in_back[di].ptr, size_out_back, mpi_type<T>::t, prev, 2*di+1, comm, &h->req[4*di+3]);
#endif
                if(gc)
                        MPI_Waitall(4, &h->req[4*di], MPI_STATUSES_IGNORE);
#elif defined(SHARK_GPI_COMM)
                assert(!pd && "Periodic domains not yet supported in GASPI!");
                int prev = domain().shiftd(di, -1, pd);
                int next = domain().shiftd(di,  1, pd);

#ifdef SHARK_STPACK
                GlobalArray<ndim, T> &ga = const_cast<GlobalArray<ndim, T> &>(*this);
                gaspi_write_notify_block(ga, ghost_out_back(di), ga, ghost_out_back(di), prev, di*2+1, notify_value);
                gaspi_write_notify_block(ga, ghost_out_front(di), ga, ghost_out_front(di), next, di*2+2, notify_value);
                if (next >= 0) h->add(di*2+1);
                if (prev >= 0) h->add(di*2+2);
#else
                const gaspi_segment_id_t seg = impl->seg;

                if (next >= 0) {
                    const gaspi_offset_t off_out_front = (g_out_front[di].ptr - ptr) * sizeof(T);
                    const gaspi_offset_t off_in_back   = (g_in_back[di].ptr - ptr) * sizeof(T);
                    const gaspi_size_t size_out_front  = g_out_front[di].r.count() * sizeof(T);
                    log_out() << "front side out:  sending " << g_out_front[di].r << std::endl;
                    log_out() << "front side out:  sending " << &g_out_front[di] << std::endl;
                    assert(g_out_front[di].full);
                    g_out_front[di].full = false;
                    SUCCESS_OR_DIE(gaspi_write_notify(seg, off_out_front, next, seg, off_in_back, size_out_front, di*2+1, notify_value, 0, GASPI_BLOCK));
                    h->add(di*2+2);
                }

                if (prev >= 0) {
                    const gaspi_offset_t off_out_back  = (g_out_back[di].ptr - ptr) * sizeof(T);
                    const gaspi_offset_t off_in_front  = (g_in_front[di].ptr - ptr) * sizeof(T);
                    const gaspi_size_t size_out_back   = g_out_back[di].r.count() * sizeof(T);
                    log_out() << "back side out: sending" << g_out_back[di].r << std::endl;
                    log_out() << "back side out:  sending " << &g_out_back[di] << std::endl;
                    assert(g_out_back[di].full);
                    g_out_back[di].full = false;
                    SUCCESS_OR_DIE(gaspi_write_notify(seg, off_out_back, prev, seg, off_in_front, size_out_back, di*2+2, notify_value, 0, GASPI_BLOCK));
                    h->add(di*2+1);
                }
#endif
               
#elif defined(SHARK_NO_COMM)

                // note: potential candidate for for_each_range
                if(pd) {
                    {
                        coords<ndim> off = ghost_back[di].adj(di, 1).lower - ghost_front[di].lower;
                        ghost_front[di].for_each([&acc,&off](coords<ndim> ii) {
                                acc(ii) = acc(ii+off);
                                });
                    }
                    {
                        coords<ndim> off = ghost_front[di].adj(di, -1).lower - ghost_back[di].lower;
                        ghost_back[di].for_each([&acc,&off](coords<ndim> ii) {
                                acc(ii) = acc(ii+off);
                                });
                    }
                }
#else
#error "No comm update"
#endif

	}

#if defined(SHARK_MPI_COMM)
	if(!gc)
		return Future<void>(std::move(h));
	else
		return Future<void>(make_unique<DoneHandle>());
#elif defined(SHARK_NO_COMM)
	return Future<void>(make_unique<DoneHandle>());
#elif defined(SHARK_GPI_COMM)
        return Future<void>(std::move(h));
#else
#error "No comm update"
#endif
}

//-- get ---

template<int ndim, typename T>
void GlobalArray<ndim,T>::get(coords_range<ndim> range, T* buf) const {
	get(range, essential_lead<ndim>(range.stride()), buf);
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::get(coords_range<ndim> range, array<size_t,ndim-1> ld, T* buf) const {
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
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
void shark::ndim::get(GlobalArray<ndim, T> &dst, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_dst) {
    assert(range_dst.counts() == range_src.counts());
#if defined(SHARK_MPI_COMM)
    typename Domain<ndim>::ProcessOverlap src_overlap(src.domain(), range_src);
    src_overlap.visit([&dst,&range_dst,&src,&range_src](int id, coords_range<ndim> i) {
#ifndef NDEBUG
            if(log_mask[verbose_rma])
               dst.log_out() << "get " << id << ": " << i << endl;
#endif
            // Target
            coords_range<ndim> t = src.domain().local(id);
            coords<ndim+1> tgtld = t.stride(src.gw);
            coord tgtoff = (i.lower - t.lower + src.gw).offset(tgtld);
            mpi_type_block<ndim,T> tgtt(i.counts(), tgtld);
     
            // Origin
            coords<ndim> o_off = i.lower - range_src.lower;  
            coords_range<ndim> local_dst = range_dst;
            local_dst.lower += o_off;
            local_dst.upper = local_dst.lower + i.counts();
            local_dst = local_dst.overlap(dst.domain().local());
            if (local_dst.count() <= 0) return;
            T *ogn = dst.pa(local_dst.lower);
            mpi_type_block<ndim,T> ognt(i.counts(), dst.ld);

            // RMA
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, id, 0, src.impl->win);
            MPI_Get(ogn, 1, ognt.t, id, tgtoff, 1, tgtt.t, src.impl->win);
            MPI_Win_unlock(id, src.impl->win);
    });
#elif defined(SHARK_GPI_COMM)
    src.log_out() << "get range_src:" << range_src << "; range_dst: " << range_dst << endl;
    // we only copy to what is local
    coords_range<ndim> local_dst = dst.domain().local().overlap(range_dst);
    // limit range_src to what is covered by local_dst
    range_src.lower += local_dst.lower - range_dst.lower;
    range_src.upper = range_src.lower + local_dst.counts();
    coords_range<ndim> off = range_src - range_dst;

    if (range_src.count() <= 0) return; //nothing to do

    typename Domain<ndim>::ProcessOverlap src_overlap(src.domain(), range_src);
    src_overlap.visit([&src,&dst,off](int id, coords_range<ndim> from) {
            //offset 'from' with same amount as 'to'
            coords_range<ndim> to = from + off;
            gaspi_read_block(src, from, dst, to, id);
            });
#elif defined(SHARK_NO_COMM)
#ifndef NDEBUG
    if(log_mask[verbose_rma])
        src.log_out() << "get " << range_src << endl;
#endif
    Access<ndim,T> dst_acc(dst);
    Access<ndim,T> src_acc(src);
    range_dst.for_each([&dst_acc,&src_acc,&range_src,&range_dst](coords<ndim> i) {
            dst_acc(i) = src_acc(i - range_dst.lower + range_src.lower);
            });
#else
#error "No comm get"
#endif
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::get(coords_range<ndim> range_src, GlobalArray<ndim, T> &dst, coords_range<ndim> range_dst) {
    auto &src = *this;
    shark::ndim::get(dst, range_src, src, range_dst);
}

//-- put ---

template<int ndim, typename T>
void shark::ndim::put(GlobalArray<ndim,T> &dst, coords_range<ndim> range_dst, GlobalArray<ndim, T> &src, coords_range<ndim> range_src) {
    assert(range_dst.counts() == range_src.counts());
#if defined(SHARK_MPI_COMM)
    typename Domain<ndim>::ProcessOverlap dst_overlap(dst.domain(), range_dst);
    dst_overlap.visit([&dst,&range_dst,&src,&range_src](int id, coords_range<ndim> i) {
            // Target
            coords_range<ndim> t = dst.domain().local(id);
            coords<ndim+1> dstld = t.stride(dst.gw);
            coord dstoff = (i.lower - t.lower + dst.gw).offset(dstld);
            mpi_type_block<ndim,T> dstt(i.counts(), dstld);

            // Origin
            coords<ndim> t_off = i.lower - range_dst.lower;  
            coords_range<ndim> local_src = range_src;
            local_src.lower += t_off;
            local_src.upper = local_src.lower + i.counts();
            local_src = local_src.overlap(src.domain().local());
            if (local_src.count() <= 0) return;
            const T *ogn = const_cast<T*>(src.pa(local_src.lower - src.domain().local().lower));
            mpi_type_block<ndim,T> ognt(i.counts(), src.ld);

            // RMA
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, id, 0, dst.impl->win);
            MPI_Put(ogn, 1, ognt.t, id, dstoff, 1, dstt.t, dst.impl->win);
            MPI_Win_unlock(id, dst.impl->win);
            });
#elif defined(SHARK_GPI_COMM)
    dst.log_out() << "put range_src:" << range_src << "; range_dst: " << range_dst << endl;
    // we only copy from what is local
    coords_range<ndim> local_src = src.domain().local().overlap(range_src);
    coords<ndim> src_offset = local_src.lower - range_src.lower;
    // limit range_dst to what is covered by local_src
    range_dst.lower += src_offset;

    if (range_dst.count() <= 0) return;

    typename Domain<ndim>::ProcessOverlap dst_overlap(dst.domain(), range_dst);
    dst_overlap.visit([&dst,&src,&range_dst,&local_src](int id, coords_range<ndim> to) {
            //offset 'from' with same amount as 'to'
            coords_range<ndim> from = local_src + to - range_dst;
            dst.log_out() << "  put from:" << from << " at rank " << dst.domain().group.procid 
            << "; to: " << to << " at rank " << id << endl;
            gaspi_write_block(src, from, dst, to, id);
            });
#elif defined(SHARK_NO_COMM)
#ifndef NDEBUG
    if(log_mask[verbose_rma])
        src.log_out() << "put " << range_src << endl;
#endif
    Access<ndim,T> dst_acc(dst);
    Access<ndim,T> src_acc(src);
    range_dst.for_each([&dst_acc,&src_acc,&range_src,&range_dst](coords<ndim> i) {
            dst_acc(i) = src_acc(i - range_dst.lower + range_src.lower);
            });
#else
#error "No comm put"
#endif
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::put(coords_range<ndim> range_dst, GlobalArray<ndim, T> &src, coords_range<ndim> range_src) {
    auto &dst = *this;
    shark::ndim::put(dst, range_dst, src, range_src);
} 

template<int ndim, typename T>
void GlobalArray<ndim,T>::put(coords_range<ndim> range, const T* buf) {
	put(range, essential_lead<ndim>(range.stride()), buf);
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::put(coords_range<ndim> range, array<size_t,ndim-1> ld, const T* buf) {
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
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

//-- accumulate ---

template<int ndim, typename T>
void shark::ndim::accumulate(GlobalArray<ndim,T> &tgt, coords_range<ndim> range_dst, GlobalArray<ndim, T> &src, coords_range<ndim> range_src) {
    assert(range_dst.counts() == range_src.counts());
#if defined(SHARK_MPI_COMM) || defined (SHARK_GPI_COMM)
    typename Domain<ndim>::ProcessOverlap tgt_overlap(tgt.domain(), range_dst);
    tgt_overlap.visit([&tgt,&range_dst,&src,&range_src](int id, coords_range<ndim> i) {
            // Target
            coords_range<ndim> t = tgt.domain().local(id);
            coords<ndim+1> tgtld = t.stride(tgt.gw);
            coord tgtoff = (i.lower - t.lower + tgt.gw).offset(tgtld);
            mpi_type_block<ndim,T> tgtt(i.counts(), tgtld);

            // Origin
            coords<ndim> t_off = i.lower - range_dst.lower;  
            coords_range<ndim> local_src = range_src;
            local_src.lower += t_off;
            local_src.upper = local_src.lower + i.counts();
            local_src = local_src.overlap(src.domain().local());
            if (local_src.count() <= 0) return;
            const T *ogn = const_cast<T*>(src.pa(local_src.lower - src.domain().local().lower));

            mpi_type_block<ndim,T> ognt(i.counts(), src.ld);

            // RMA
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, id, 0, tgt.impl->win);
            MPI_Accumulate(ogn, 1, ognt.t, id, tgtoff, 1, tgtt.t, MPI_SUM, tgt.impl->win);
            MPI_Win_unlock(id, tgt.impl->win);
            });
#elif defined(SHARK_NO_COMM)
#ifndef NDEBUG
    if(log_mask[verbose_rma])
        src.log_out() << "accumulate " << range_src << endl;
#endif
    Access<ndim,T> tgt_acc(tgt);
    Access<ndim,T> src_acc(src);
    range_dst.for_each([&tgt_acc,&src_acc,&range_src](coords<ndim> i) {
            tgt_acc(i) += src_acc(i - range_src.lower);
            });
#else
#error "No comm accumulate"
#endif
}

template<int ndim, typename T>
void GlobalArray<ndim,T>::accumulate(coords_range<ndim> range_dst, GlobalArray<ndim, T> &src, coords_range<ndim> range_src) {
    shark::ndim::accumulate(*this, range_dst, src, range_src);
} 

template<int ndim,typename T> template<typename>
void GlobalArray<ndim,T>::dump(coords_range<ndim> range, std::string filename) {
#ifndef NDEBUG
        if(log_mask[verbose_rma])
                this->log_out() << "dump to " << filename << endl;
#endif
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

        // type for part to read from
        coords_range<ndim> src = domain().local().overlap(range);
        coords<ndim+1> srcld = src.stride(gw);
        const T* buf = pa(src.lower);
        mpi_type_block<ndim,T> srct(src.counts(), srcld);
        
        // where to write to
        coords_range<ndim> dst = region().overlap(range);
        coords<ndim+1> dstld = dst.stride(); // no ghost regions in file
        coord dstoff = dst.lower.offset(dstld); // offset wrt to global array
        
        MPI_File_write_at(fh, dstoff, const_cast<T*>(buf), 1, srct.t, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);
#elif defined(SHARK_NO_COMM)
        FILE *fd = fopen(filename.c_str(), "w");
        Access<ndim,T> acc(*this);
	range.for_each([&acc,fd](coords<ndim> ii) {
                fwrite(&acc(ii), sizeof(T), 1, fd);
	});
        fclose(fd);
#else
#error "No comm dump"
#endif
}

template<int ndim,typename T>
GADest<ndim,T> GlobalArray<ndim,T>::region(coords_range<ndim> r, bool pack) {
	return GADest<ndim,T>(*this, r, pack);
}

template<int ndim,typename T>
GADest<ndim,T>::GADest(GlobalArray<ndim,T>& ga, coords_range<ndim> r, bool pack)
    : ga(ga), region(r), pack(pack) { }

template<int ndim,typename T>
GADest<ndim,T>::~GADest() { }

template<int ndim,typename T>
GADest<ndim,T>::GADest(const GADest& gad): ga(gad.ga), region(gad.region), pack(gad.pack) { }

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

#define SYMBDT(d,T) template void  shark::ndim::GlobalArray<d,T>::dump(coords_range<d> range, std::string);
#include "inst_dimtype"
#undef SYMBDT

#define SYMBDT(d,T) template ostream& shark::ndim::operator<< <d,T>(ostream& os, const shark::ndim::GlobalArray<d,T> &ga);
#include "inst_dimtype"
#undef SYMBDT

#undef SYMB_ARITH

#define SYMBDT(d,T) template class shark::ndim::GADest<d,T>; 
#include "inst_dimtype"
#undef SYMBDT
	
#define SYMBDT(d,T) template class shark::ndim::GARef<d,T>; 
#include "inst_dimtype"
#undef SYMBDT	

#define SYMBDT(d,T) template class shark::ndim::GABuf<d,T>; 
#include "inst_dimtype"
#undef SYMBDT
