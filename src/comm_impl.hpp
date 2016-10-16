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
#include <shark/group.hpp>
#include <shark/globals.hpp>
#include <shark/common.hpp>

#include <set>
#include <valarray>
#include <complex>
#include <cassert>

#if defined(SHARK_MPI_COMM) || defined(SHARK_GPI_COMM)

#define OMPI_SKIP_MPICXX              // OpenMPI: only C bindings
#define MPICH_SKIP_MPICXX             // MPICH: only C bindings
#define MPI_NO_CPPBIND                // MPT: only C bindings
#include <mpi.h>

#ifdef SHARK_GPI_COMM
#include <GASPI.h>
#include <GASPI_Ext.h>
#ifdef NDEBUG
#define success_or_die(f...) f
#else
#define success_or_die(f...) assert (f == GASPI_SUCCESS); 
#endif
#endif

namespace shark {
	class GroupImpl {
	public:
		MPI_Comm comm;

#ifdef SHARK_GPI_COMM
                gaspi_group_t grp;
#endif
		int rank();
		int size();

#ifdef SHARK_GPI_COMM
                ~GroupImpl();

                static std::set<unsigned> segs_used;
                static void free_seg_id(int i) {
                    assert(segs_used.find(i) != segs_used.end());
                    segs_used.erase(i);
                }
                static int next_seg_id() {
                    gaspi_number_t segment_max;
                    gaspi_segment_max(&segment_max);
                    for(unsigned i=0; i<segment_max; ++i) {
                        if (segs_used.find(i) == segs_used.end()) {
                            segs_used.insert(i);
                            return i;
                        }
                    }
                    assert(false);
                    return -1;
                }
#endif

	};

	namespace ndim {
		template<int ndim, typename T>
		class GAImpl {
		public:
			MPI_Win win;
			MPI_Datatype ghost[ndim];
#ifdef SHARK_GPI_COMM
			gaspi_segment_id_t seg;
#endif
		};
		template<int ndim>
		class SAImpl {
		public:
#ifdef SHARK_GPI_COMM
			gaspi_segment_id_t seg;
#endif 
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

#if defined(SHARK_GPI_COMM)

#define SUCCESS_OR_DIE(f...) \
do  { \
  gaspi_number_t queue_size, queue_max; \
  gaspi_queue_size_max(&queue_max); \
  gaspi_queue_size(0, &queue_size); \
  if(queue_size > queue_max - 1) gaspi_wait(0, GASPI_BLOCK); \
  \
  const gaspi_return_t r = f; \
  if (r != GASPI_SUCCESS) { \
    gaspi_printf ("Error: '%s' [%s:%i]: %i\n", #f, __FILE__, __LINE__, r); \
    abort(); \
  } \
} while (0);


namespace shark {

	namespace ndim {

                template<int ndim, typename A, typename B>
                void gaspi_write_block(A &local, coords_range<ndim> range_local, B &remote, coords_range<ndim> range_remote, int rank) {
                    SHARK_COUNTER("gaspi_write_block");
                    if (rank == MPI_PROC_NULL) return;
                    if (range_local.count() <= 0) return;
                    if (range_remote.count() <= 0) return;

                    local.log_out() << "Writing block to rank " << rank << "; from " << range_local << "; to " << range_remote 
                        << "; counts " << range_local.counts() << std::endl;

                    assert (remote.domain().local(rank, remote.ghost_width()).contains(range_remote) && 
                            "gaspi_write_block: remote domain does not contain remote range");

                    assert (local.domain().local(local.ghost_width()).contains(range_local) &&
                            "gaspi_write_block: local domain does not contain local range");

                    coords<ndim> offset = range_remote.lower - range_local.lower;

                    range_local.for_block([&local, &remote, offset, rank](coords<ndim> i, coord count) {
                            gaspi_segment_id_t  seg_local = local.impl->seg;
                            gaspi_segment_id_t  seg_remote = remote.impl->seg;

                            coords<ndim> from = i - local.domain().local().lower;
                            coords<ndim> to   = i + offset - remote.domain().local(rank).lower;

#ifndef NDEBUG
                            local.log_out() << "  Writing to rank " << rank << "; from " << from << "; to " << to << "; count " << count << std::endl;
                            local.log_out() << "    Byte offset from: " << local.byte_offset(from)  << std::endl;
                            local.log_out() << "    Byte offset to: " << remote.byte_offset(to)  << std::endl;
                            local.log_out() << "    Size: " << count * sizeof(typename A::value_type) << std::endl;
#endif

                            {SHARK_COUNTER("gaspi_write"); SUCCESS_OR_DIE(gaspi_write(seg_local, local.byte_offset(from), rank, seg_remote, remote.byte_offset(to), count * sizeof(typename A::value_type), 0, GASPI_BLOCK));
                            }
                            });

                    local.log_out() << "Done..." << std::endl;
                }

                template<int ndim, typename A, typename B>
                void gaspi_read_block(A &local, coords_range<ndim> range_local, B &remote, coords_range<ndim> range_remote, int rank) {
                    if (rank == MPI_PROC_NULL) return;
                    if (range_local.count() <= 0) return;
                    if (range_remote.count() <= 0) return;
                    assert(range_local.counts() == range_remote.counts());
                    
                    local.log_out() << "Reading block from rank " << rank << "; to " << range_local << "; from " << range_remote 
                                  << "; counts " << range_local.counts() << std::endl;

                    coords<ndim> offset = range_local.lower - range_remote.lower;

                    range_local.for_block([&local, &remote, offset, rank](coords<ndim> i, coord count) {
                            gaspi_segment_id_t  seg_local = local.impl->seg;
                            gaspi_segment_id_t  seg_remote = remote.impl->seg;

                            coords<ndim> to = i - local.domain().local().lower;
                            coords<ndim> from   = i + offset - remote.domain().local(rank).lower;

#ifndef NDEBUG
                            local.log_out() << "  Reading from rank " << rank << "; to " << to << "; from " << from << "; count " << count << std::endl;
#endif

                            SUCCESS_OR_DIE(gaspi_read(seg_local, local.byte_offset(to), rank, seg_remote, remote.byte_offset(from), count * sizeof(typename A::value_type), 0, GASPI_BLOCK));
                    });

                     local.log_out() << "Done..." << std::endl;
                }

                template<int ndim, typename A, typename B>
                void gaspi_write_notify_block(A &local, coords_range<ndim> range_local, B &remote, coords_range<ndim> range_remote, int rank, 
                        int notification_id, int notification_value = 1) {
                    if (rank == MPI_PROC_NULL) return;

                    gaspi_write_block(local, range_local, remote, range_remote, rank);
                    SUCCESS_OR_DIE(gaspi_notify(remote.impl->seg, rank, notification_id, notification_value, 0, GASPI_BLOCK));
                }
        }
}
#endif

#elif defined(SHARK_NO_COMM)

#include <cstddef>   // std::size_t

namespace shark {
	void* mem_alloc(std::size_t size);
	void mem_free(void* p);

	class GroupImpl { };

	namespace ndim {
		template<int ndim, typename T>
		class GAImpl {
                public:
                    GAImpl() {} 
                };
		template<int ndim>
		class SAImpl { };
	}
}

#else

#error "No comm defined"

#endif

#endif
