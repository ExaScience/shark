#ifndef __SHARK_MPI_IMPL_HPP
#define __SHARK_MPI_IMPL_HPP

#include <shark/vec.hpp>

#define OMPI_SKIP_MPICXX              // OpenMPI: only C bindings
#define MPICH_SKIP_MPICXX             // MPICH: only C bindings
#define MPI_NO_CPPBIND                // MPT: only C bindings
#include <mpi.h>

namespace shark {
	class GroupImpl {
	public:
		MPI_Comm comm;
		int rank();
		int size();
	};

	namespace ndim {
		template<int ndim>
		class GAImpl {
		public:
			MPI_Win win;
			MPI_Datatype ghost[ndim];
		};
	}

	template<typename T>
	struct mpi_type {
		static const MPI_Datatype t;
		INLINE static int count(const T& proto);
		INLINE static void* address(T& object);
	};

	// Inline function implementations

	template<typename T>
	inline int mpi_type<T>::count(const T&) {
		return 1;
	}

	template<typename T>
	inline void* mpi_type<T>::address(T& object) {
		return &object;
	}

	template<int ndim,typename T>
	struct mpi_type<shark::ndim::vec<ndim,T>> {
		static MPI_Datatype t;
		static void init();
		static void destroy();
		INLINE static int count(const shark::ndim::vec<ndim,T>& proto);
		INLINE static void* address(shark::ndim::vec<ndim,T>& object);
	};

	template<int ndim,typename T>
	inline int mpi_type<shark::ndim::vec<ndim,T>>::count(const shark::ndim::vec<ndim,T>&) {
		return 1;
	}

	template<int ndim,typename T>
	inline void* mpi_type<shark::ndim::vec<ndim,T>>::address(shark::ndim::vec<ndim,T>& object) {
		return &object.val;
	}
}

#define SHARK_MPI_INT_INST \
	SYMB(int) \
	SYMB(long)
#define SHARK_MPI_FP_INST \
	SYMB(float) \
	SYMB(double) \
	SYMB(long double) 

#if (MPI_VERSION > 2 || MPI_VERSION == 2 && MPI_SUBVERSION >= 2)
#define SHARK_MPI_COMP_INST \
	SYMB(std::complex<float>) \
	SYMB(std::complex<double>) \
	SYMB(std::complex<long double>)
#else
#define SHARK_MPI_COMP_INST
#endif

#endif
