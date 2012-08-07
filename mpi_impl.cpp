
#include "mpi_impl.hpp"

using namespace shark;
using namespace shark::ndim;

template<> const MPI_Datatype mpi_type<int>::t = MPI_INT;
template<> const MPI_Datatype mpi_type<long>::t = MPI_LONG;
template<> const MPI_Datatype mpi_type<float>::t = MPI_FLOAT;
template<> const MPI_Datatype mpi_type<double>::t = MPI_DOUBLE;

template<int ndim,typename T> MPI_Datatype mpi_type<vec<ndim,T>>::t;

template<int ndim,typename T>
void mpi_type<vec<ndim,T>>::init() {
	MPI_Type_contiguous(ndim, mpi_type<T>::t, &t);
	MPI_Type_commit(&t);
}

template<int ndim,typename T>
void mpi_type<vec<ndim,T>>::destroy() {
	MPI_Type_free(&t);
}

// Set-up instantiations

#define SYMB(d,T) template struct mpi_type<vec<d,T>>; 
#include "inst_dimtype"
#undef SYMB


