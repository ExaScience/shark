/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <utility>                     // std::move
#include <iostream>

#include <shark/group.hpp>
#include "comm_impl.hpp"

using namespace std;
using namespace shark;

#if defined(SHARK_MPI_COMM) || defined (SHARK_GPI_COMM)

#ifdef SHARK_GPI_COMM
std::set<unsigned> GroupImpl::segs_used;

GroupImpl::~GroupImpl() {
    assert(segs_used.empty());
}
#endif

int GroupImpl::rank() {
	int rank;
	MPI_Comm_rank(comm, &rank);
	return rank;
}

int GroupImpl::size() {
	int size;
	MPI_Comm_size(comm, &size);
	return size;
}

#endif

unique_ptr<Group> Group::w;
unique_ptr<Group> Group::s;

#if defined(SHARK_MPI_COMM) || defined (SHARK_GPI_COMM)
Group::Group(unique_ptr<GroupImpl>&& oimpl): impl(move(oimpl)), procid(impl->rank()), nprocs(impl->size()) { }
#elif defined(SHARK_NO_COMM)
Group::Group(unique_ptr<GroupImpl>&& oimpl): impl(move(oimpl)), procid(0), nprocs(1) { }
#endif

Group::~Group() { }

bool Group::operator==(const Group& other) const {
	return this == &other;
}

void Group::sync() const {
#if defined(SHARK_MPI_COMM) || defined (SHARK_GPI_COMM)
	MPI_Barrier(impl->comm);
#if  defined (SHARK_GPI_COMM)
        gaspi_wait(0, GASPI_BLOCK);
        gaspi_barrier(impl->grp,GASPI_BLOCK);
#endif
#endif
}

template<typename T>
T Group::external_sum(T&& sum) const
{
#if defined(SHARK_MPI_COMM) || defined (SHARK_GPI_COMM)
	MPI_Allreduce(MPI_IN_PLACE, mpi_type<T>::address(sum), mpi_type<T>::count(sum), mpi_type<T>::t, MPI_SUM, impl->comm);
	return sum;
#elif defined(SHARK_NO_COMM)
	return sum;
#else
#error "No comm external_sum"
#endif
}

#if defined(SHARK_MPI_COMM) || defined (SHARK_GPI_COMM)
namespace {
	class SumHandle: public Handle {
		friend class shark::Group;
		MPI_Request req;
	public:
		SumHandle();
		virtual ~SumHandle();
		virtual bool test();
		virtual void wait();
	};
}

SumHandle::SumHandle() {
}

SumHandle::~SumHandle() {
}

bool SumHandle::test() {
	int flag;
	MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
	return flag;
}

void SumHandle::wait() {
	MPI_Wait(&req, MPI_STATUS_IGNORE);
}
#endif

template<typename T>
Future<T> Group::external_isum(T&& sum) const
{
#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)
	auto val = make_unique<T>(std::move(sum));
	auto h = make_unique<SumHandle>();
	MPI_Iallreduce(MPI_IN_PLACE, mpi_type<T>::address(*val), mpi_type<T>::count(*val), mpi_type<T>::t, MPI_SUM, impl->comm, &h->req);
	return Future<T>(std::move(h), std::move(val));
#elif defined(SHARK_NO_COMM)
	return Future<T>(make_unique<DoneHandle>(), make_unique<T>(std::move(sum)));
#else
#error "No comm external_isum"
#endif
}

// Explicit instantiations

#define SYMBT(T) template T Group::external_sum(T&&) const;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#undef SYMBT

#define SYMBT(T) template Future<T> Group::external_isum(T&&) const;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#undef SYMBT

