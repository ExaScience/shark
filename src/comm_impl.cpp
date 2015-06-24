/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include "comm_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

#if defined(SHARK_MPI_COMM)

template<> const MPI_Datatype mpi_type<char>::t = MPI_CHAR;
template<> const MPI_Datatype mpi_type<int>::t = MPI_INT;
template<> const MPI_Datatype mpi_type<long>::t = MPI_LONG;
template<> const MPI_Datatype mpi_type<float>::t = MPI_FLOAT;
template<> const MPI_Datatype mpi_type<double>::t = MPI_DOUBLE;
template<> const MPI_Datatype mpi_type<long double>::t = MPI_LONG_DOUBLE;
template<> const MPI_Datatype mpi_type<std::complex<float>>::t = MPI_C_COMPLEX;
template<> const MPI_Datatype mpi_type<std::complex<double>>::t = MPI_C_DOUBLE_COMPLEX;
template<> const MPI_Datatype mpi_type<std::complex<long double>>::t = MPI_C_LONG_DOUBLE_COMPLEX;

template<int ndim> const MPI_Datatype mpi_type<coords_range<ndim>>::t = mpi_type<coord>::t;

template<int ndim,typename T> const MPI_Datatype& mpi_type<vec<ndim,T>>::t(mpi_type<T>::t);
MPI_Datatype mpi_type<part_position>::t;
template<int ndim> MPI_Datatype mpi_type<part<ndim>>::t;

template<typename T> const MPI_Datatype& mpi_type<valarray<T>>::t(mpi_type<T>::t);

const MPI_Datatype mpi_type<test_result>::t = MPI_LONG;

template<int ndim, typename T>
MPI_Datatype mpi_type_block<ndim,T>::type(coords<ndim> n, array<size_t,ndim-1> eld) {
	MPI_Datatype t;
	MPI_Aint lb, extent;
	if(mpi_type<T>::count() != 1)
		MPI_Type_contiguous(mpi_type<T>::count(), mpi_type<T>::t, &t);
	else
		MPI_Type_dup(mpi_type<T>::t, &t);
	MPI_Type_get_extent(t, &lb, &extent);
	for(int di = ndim-1; di >= 0; di--) {
		MPI_Datatype tmp;
		tmp = t;
		MPI_Type_create_hvector(n[di], 1, di == ndim-1 ? extent : eld[di]*extent, tmp, &t);
		MPI_Type_free(&tmp);
	}
	MPI_Type_commit(&t);
	return t;
}

template<int ndim, typename T>
mpi_type_block<ndim,T>::mpi_type_block(coords<ndim> n, array<size_t,ndim-1> eld): t(type(n,eld)) {
}

template<int ndim, typename T>
mpi_type_block<ndim,T>::mpi_type_block(coords<ndim> n, coords<ndim+1> ld): t(type(n,essential_lead<ndim>(ld))) {
}

template<int ndim, typename T>
mpi_type_block<ndim,T>::~mpi_type_block() {
	MPI_Datatype tmp = t;
	MPI_Type_free(&tmp);
}

// Set-up instantiations

#define SYMBT(T) template struct shark::mpi_type<T>;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#include "comm_other_inst"
#undef SYMBT

#define SYMBDT(d,T) template class shark::ndim::mpi_type_block<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

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
