
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

#include "types"

#define SYMB(T) template T Group::external_sum(T&&) const;
SHARK_COMM_INT_INST
SHARK_COMM_FP_INST
SHARK_COMM_CPLX_INST
#undef SYMB

#define SYMB(T) template valarray<T> Group::external_sum(valarray<T>&&) const;
SHARK_COMM_INT_INST
SHARK_COMM_FP_INST
SHARK_COMM_CPLX_INST
#undef SYMB

#define SYMB(d,T) template ndim::vec<d,T> Group::external_sum(ndim::vec<d,T>&&) const; 
#include "inst_dimtype"
#undef SYMB

#ifdef ECL_ASYNC
#define SYMB(T) template T Group::external_isum(T&&) const;
SHARK_COMM_INT_INST
SHARK_COMM_FP_INST
SHARK_COMM_CPLX_INST
#undef SYMB
#endif

