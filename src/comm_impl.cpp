/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include "comm_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

#if defined(SHARK_MPI_COMM)

template<> const MPI_Datatype mpi_type<int>::t = MPI_INT;
template<> const MPI_Datatype mpi_type<long>::t = MPI_LONG;
template<> const MPI_Datatype mpi_type<float>::t = MPI_FLOAT;
template<> const MPI_Datatype mpi_type<double>::t = MPI_DOUBLE;
template<> const MPI_Datatype mpi_type<long double>::t = MPI_LONG_DOUBLE;
template<> const MPI_Datatype mpi_type<std::complex<float>>::t = MPI_C_COMPLEX;
template<> const MPI_Datatype mpi_type<std::complex<double>>::t = MPI_C_DOUBLE_COMPLEX;
template<> const MPI_Datatype mpi_type<std::complex<long double>>::t = MPI_C_LONG_DOUBLE_COMPLEX;

template<int ndim,typename T> const MPI_Datatype mpi_type<vec<ndim,T>>::t = mpi_type<T>::t;
template<int ndim,typename T> const MPI_Datatype mpi_type<part<ndim,T>>::t = mpi_type<T>::t;

template<typename T> const MPI_Datatype mpi_type<valarray<T>>::t = mpi_type<T>::t;

const MPI_Datatype mpi_type<test_result>::t = MPI_LONG;

// Set-up instantiations

#define SYMBT(T) template struct shark::mpi_type<T>;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#include "comm_other_inst"
#undef SYMBT

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
