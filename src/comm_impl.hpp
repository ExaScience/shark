/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_COMM_IMPL_HPP
#define __SHARK_COMM_IMPL_HPP

#include <shark/coords_range.hpp>
#include <shark/vec.hpp>
#include <shark/part.hpp>
#include <valarray>
#include <complex>

#if defined(SHARK_MPI_COMM)

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

	/*
	 * mpi_type for predefined MPI types
	 */
	template<typename T>
	struct mpi_type {
		static const MPI_Datatype t;
		INLINE static int count();
		INLINE static int count(const T& proto);
		INLINE static void* address(T& object);
	};

	// Inline function implementations
	
	template<typename T>
	inline int mpi_type<T>::count() {
		return 1;
	}

	template<typename T>
	inline int mpi_type<T>::count(const T&) {
		return 1;
	}

	template<typename T>
	inline void* mpi_type<T>::address(T& object) {
		return &object;
	}

	/*
	 * mpi_type for coord_range
	 */
	template<int ndim>
	struct mpi_type<shark::ndim::coords_range<ndim>> {
		static const MPI_Datatype t;
		INLINE static int count();
	};

	template<int ndim>
	inline int mpi_type<shark::ndim::coords_range<ndim>>::count() {
		return ndim*2;
	}

	/*
	 * mpi_type for vec
	 */
	template<int ndim,typename T>
	struct mpi_type<shark::ndim::vec<ndim,T>> {
		static const MPI_Datatype& t;
		INLINE static int count();
		INLINE static int count(const shark::ndim::vec<ndim,T>& proto);
		INLINE static void* address(shark::ndim::vec<ndim,T>& object);
	};

	template<int ndim,typename T>
	inline int mpi_type<shark::ndim::vec<ndim,T>>::count() {
		return ndim;
	}

	template<int ndim,typename T>
	inline int mpi_type<shark::ndim::vec<ndim,T>>::count(const shark::ndim::vec<ndim,T>&) {
		return ndim;
	}

	template<int ndim,typename T>
	inline void* mpi_type<shark::ndim::vec<ndim,T>>::address(shark::ndim::vec<ndim,T>& object) {
		return &object.val;
	}

	/*
	 * mpi_type for part_position
	 */
	template<>
	struct mpi_type<shark::ndim::part_position> {
		static MPI_Datatype t;
		INLINE static int count();
	};

	inline int mpi_type<shark::ndim::part_position>::count() {
		return 1;
	}

	/*
	 * mpi_type for part
	 */
	template<int ndim>
	struct mpi_type<shark::ndim::part<ndim>> {
		static MPI_Datatype t;
		INLINE static int count();
	};

	template<int ndim>
	inline int mpi_type<shark::ndim::part<ndim>>::count() {
		return 1;
	}

	/*
	 * mpi_type for valarray
	 */
	template<typename T>
	struct mpi_type<std::valarray<T>> {
		static const MPI_Datatype& t;
		INLINE static int count(const std::valarray<T>& proto);
		INLINE static void* address(std::valarray<T>& object);
	};

	template<typename T>
	inline int mpi_type<std::valarray<T>>::count(const std::valarray<T>& proto) {
		return static_cast<int>(proto.size());
	}

	template<typename T>
	inline void* mpi_type<std::valarray<T>>::address(std::valarray<T>& object) {
		return &object[0];
	}

	template<>
	struct mpi_type<test_result> {
		static const MPI_Datatype t;
		INLINE static int count();
		INLINE static int count(const test_result& proto);
		INLINE static void* address(test_result& object);
	};

	inline int mpi_type<test_result>::count() {
		return 2;
	}

	inline int mpi_type<test_result>::count(const test_result&) {
		return 2;
	}

	inline void* mpi_type<test_result>::address(test_result& object) {
		return &object;
	}

	namespace ndim {

		template<int ndim, typename T>
		class mpi_type_block {
			static MPI_Datatype type(coords<ndim> n, std::array<std::size_t,ndim-1> eld);
		public:
			const MPI_Datatype t;
			mpi_type_block(coords<ndim> n, coords<ndim+1> ld);
			mpi_type_block(coords<ndim> n, std::array<std::size_t,ndim-1> eld);
			~mpi_type_block();
		};

	}
}

#elif defined(SHARK_NO_COMM)

#include <cstddef>   // std::size_t

namespace shark {
	void* mem_alloc(std::size_t size);
	void mem_free(void* p);

	class GroupImpl { };

	namespace ndim {
		template<int ndim>
		class GAImpl { };
	}
}

#else

#error "No comm defined"

#endif

#endif
