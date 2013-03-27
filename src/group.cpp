/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <utility>                     // std::move

#include <shark/group.hpp>
#include "comm_impl.hpp"

using namespace std;
using namespace shark;

#if defined(SHARK_MPI_COMM)

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

#if defined(SHARK_MPI_COMM)
Group::Group(unique_ptr<GroupImpl>&& oimpl): impl(move(oimpl)), procid(impl->rank()), nprocs(impl->size()) { }
#elif defined(SHARK_NO_COMM)
Group::Group(unique_ptr<GroupImpl>&& oimpl): impl(move(oimpl)), procid(0), nprocs(1) { }
#endif

Group::~Group() { }

bool Group::operator==(const Group& other) const {
	return this == &other;
}

void Group::sync() const {
#if defined(SHARK_MPI_COMM)
	MPI_Barrier(impl->comm);
#elif defined(SHARK_NO_COMM)
#endif
}

template<typename T>
T Group::external_sum(T&& sum) const {
#if defined(SHARK_MPI_COMM)
	MPI_Allreduce(MPI_IN_PLACE, mpi_type<T>::address(sum), mpi_type<T>::count(sum), mpi_type<T>::t, MPI_SUM, impl->comm);
#endif
	return sum;
}

#ifdef ECL_ASYNC
template<typename T>
Future<T> Group::external_isum(T&& sum) const {
	Future<T> f(false, std::move(sum));
	MPIX_Iallreduce(MPI_IN_PLACE, mpi_type<T>::address(*f.val), mpi_type<T>::count(*f.val), mpi_type<T>::t, MPI_SUM, impl->comm, &f.h->r);
	return f;
}
#endif

// Explicit instantiations

#define SYMBT(T) template T Group::external_sum(T&&) const;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#undef SYMBT

#ifdef ECL_ASYNC
#define SYMBT(T) template T Group::external_isum(T&&) const;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#undef SYMBT
#endif

