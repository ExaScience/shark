
#include <utility>                     // std::move

#include <shark/group.hpp>
#include "comm_impl.hpp"

using namespace std;
using namespace shark;

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

unique_ptr<Group> Group::w;

Group::Group(unique_ptr<GroupImpl>&& oimpl): impl(move(oimpl)), procid(impl->rank()), nprocs(impl->size()) { }

Group::~Group() { }

bool Group::operator==(const Group& other) const {
	return this == &other;
}

void Group::sync() const {
	MPI_Barrier(impl->comm);
}

template<typename T>
T Group::external_sum(T&& sum) const {
	MPI_Allreduce(MPI_IN_PLACE, mpi_type<T>::address(sum), mpi_type<T>::count(sum), mpi_type<T>::t, MPI_SUM, impl->comm);
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

#define SYMB(T) template T Group::external_sum(T&&) const;
SHARK_MPI_INT_INST
SHARK_MPI_FP_INST
SHARK_MPI_COMP_INST
#undef SYMB

#ifdef ECL_ASYNC
#define SYMB(T) template T Group::external_isum(T&&) const;
SHARK_MPI_INT_INST
SHARK_MPI_FP_INST
SHARK_MPI_COMP_INST
#undef SYMB
#endif

