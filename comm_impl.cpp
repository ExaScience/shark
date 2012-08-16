
#include "comm_impl.hpp"

#include <complex>

using namespace std;
using namespace shark;
using namespace shark::ndim;

#if defined(SHARK_MPI_COMM)

#include <complex>

template<> const MPI_Datatype mpi_type<int>::t = MPI_INT;
template<> const MPI_Datatype mpi_type<long>::t = MPI_LONG;
template<> const MPI_Datatype mpi_type<float>::t = MPI_FLOAT;
template<> const MPI_Datatype mpi_type<double>::t = MPI_DOUBLE;
template<> const MPI_Datatype mpi_type<long double>::t = MPI_LONG_DOUBLE;
#if (MPI_VERSION > 2 || MPI_VERSION == 2 && MPI_SUBVERSION >= 2)
template<> const MPI_Datatype mpi_type<std::complex<float>>::t = MPI_C_COMPLEX;
template<> const MPI_Datatype mpi_type<std::complex<double>>::t = MPI_C_DOUBLE_COMPLEX;
template<> const MPI_Datatype mpi_type<std::complex<long double>>::t = MPI_C_LONG_DOUBLE_COMPLEX;
#endif

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

#elif defined(SHARK_NO_COMM)

#include <cstdlib>

void* shark::mem_alloc(size_t size) {
	return valloc(size);
}

void shark::mem_free(void* p) {
	free(p);
}

#else

#error "No comm defined"

#endif
