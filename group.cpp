
#include <utility>                     // std::move

#include <shark/group.hpp>
#include "mpi_impl.hpp"

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


