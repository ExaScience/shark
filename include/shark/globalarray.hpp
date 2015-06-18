/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_GLOBALARRAY_HPP
#define __SHARK_GLOBALARRAY_HPP

#include <array>                       // std::array
#include <cstddef>                     // std::size_t
#include <memory>                      // std::unique_ptr
#include <cassert>                     // assert
#include "common.hpp"
#include "boundary.hpp"
#include "coords.hpp"
#include "coords_range.hpp"
#include "future.hpp"
//#include <mpi.h>

namespace shark {

	namespace ndim {

		template<int>
		class GAImpl;

		template<int,typename>
		class GADest;

		template<int,typename>
		class GARef;

		/**
		 * A global array that is partitioned across a process group according to a domain
		 */
		template<int ndim,typename T>
		class GlobalArray
		{
			friend class Access<ndim,T>;
		
		public:
			typedef GARef<ndim,T> storage;
			typedef Access<ndim,T> accessor;
			static const /*constexpr*/ int number_of_dimensions = ndim;
			typedef std::array<Boundary<ndim,T>,ndim> bounds;

		private:
			// before allocate
			const Domain<ndim>* dom;
			coords<ndim> gw;
			bool gc;
			bounds bd;

			// after allocate
			T* ptr;
			std::unique_ptr<GAImpl<ndim>> impl;
			coords<ndim+1> ld;
			std::array<coords_range<ndim>,ndim> ghost_back;
			std::array<coords_range<ndim>,ndim> ghost_front;

			// extra
			mutable int lc;

			void allocate();
			void deallocate();
			void reset();

			INLINE T& da(coords<ndim> i) const;

			std::ostream& log_out() const;

		public:
			/**
			 * The domain
			 */
			INLINE const Domain<ndim>& domain() const;

			/**
			 * The number of ghost cells for each dimension.
			 */
			INLINE coords<ndim> ghost_width() const;

			/**
			 * Whether corner ghost cells are maintained.
			 */
			INLINE bool ghost_corners() const;

			/**
			 * Returns whether this global array is active.
			 */
			INLINE /*explicit*/ operator bool() const;

			/**
			 * The total/inner/outer region of elements
			 */
			INLINE coords_range<ndim> region() const;
			INLINE coords_range<ndim> inner() const;
			INLINE coords_range<ndim> outer_front(int) const;
			INLINE coords_range<ndim> outer_back(int) const;

			/**
			 * Select a region of elements as destination
			 */
			GADest<ndim,T> region(coords_range<ndim> r);

			/**
			 * Construct a GlobalArray (collective).
			 * The global array will not be active until it is assigned.
			 */
			GlobalArray();

			/**
			 * Construct a GlobalArray (collective).
			 * @param domain the domain for the distribution of data
			 * @param ghost_width number of ghost cells for each dimension
			 * @param ghost_corners whether to maintain the corner ghost cells
			 */
			GlobalArray(const Domain<ndim>& domain, coords<ndim> ghost_width = coords<ndim>(), bool ghost_corners = false, bounds bd = bounds());

			GlobalArray(const GlobalArray<ndim,T>& other, bool copy);

			/**
			 * Destruct a GlobalArray (collective).
			 * If active, the memory of the global array is released
			 */
			~GlobalArray();

			// Move semantics
			GlobalArray(const GlobalArray<ndim,T>& other) = delete;

			GlobalArray(GlobalArray<ndim,T>&& other);
			GlobalArray<ndim,T>& operator=(GlobalArray<ndim,T>&& other);

			// Copy from source
			GlobalArray<ndim,T>& operator=(const GlobalArray<ndim,T>& other);
			template<typename S>
			GlobalArray<ndim,T>& operator=(const S&);

			/**
			 * Redistribute the GlobalArray to have a Domain.
			 * @param domain the new domain, \verbatim domain().equiv(domain) \enverbatim
			 */
			void reshape(const Domain<ndim>& domain);

			/**
			 * Update the ghost cells with data from their original locations (collective).
			 * RMA operations cannot overlap with local access.
			 *
			 * iupdate is a non-blocking variant if circumstances permit this (aync communication
			 * support, no ghost corners)
			 *
			 * @param k update phase (only used for general boundaries)
			 */
			void update(long k=0) const;
			Future<void> iupdate(long k=0) const;

			/**
			 * Get remote range (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to retrieve
			 * @param ld the strides to use for buf (default: determined by range)
			 * @param buf the target buffer
			 */

			void get(coords_range<ndim> range, T* buf) const;
			void get(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, T* buf) const;

			/**
			 * Put remote range (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to retrieve
			 * @param ld the strides to use for buf (default: determined by range)
			 * @param buf the target buffer
			 */
			void put(coords_range<ndim> range, const T* buf);
			void put(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, const T* buf);

			/**
			 * Increase remote range with local values (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to update
			 * @param ld the strides to use for buf (default: determined by range)
			 * @param buf the source buffer
			 */
			template<typename = void>
			void accumulate(coords_range<ndim> range, const T* buf);
			template<typename = void>
			void accumulate(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, const T* buf);

			/**
			 * Selective gather of remote values (collective). This fills all available positions in a
			 * local sparse array with their values from this global array.
			 * RMA operations cannot overlap with local access.
			 * @param sa the sparse array to fill
			 *
			 * igather is a non-blocking variant if circumstances permit this (aync communication
			 * support)
			 */
			template<typename = void>
			void gather(SparseArray<ndim,T>& sa) const;
			template<typename = void>
			Future<void> igather(SparseArray<ndim,T>& sa) const;

			/**
			 * Selective scatter to remote values (collective). This sends out all available values of a
			 * local sparse array to the respective positions in this global array. The values are accumulated
			 * with the original value at the destination.
			 * RMA operations cannot overlap with local access.
			 * @param sa the sparse array whose values to send out
			 *
			 * iscatterAcc is a non-blocking variant if circumstances permit this (aync communication
			 * support)
			 */
			template<typename = void>
			void scatterAcc(const SparseArray<ndim,T>& sa);
			template<typename = void>
			Future<void> iscatterAcc(const SparseArray<ndim,T>& sa);

			/**
			 * Dump values of global array to a file
			 */			template<typename = void>
			void dump(coords_range<ndim> range, std::string filename);
		};

		template<int ndim, typename T>
		class GADest {
			friend class GlobalArray<ndim,T>;
			GlobalArray<ndim,T>& ga;
			coords_range<ndim> r;
			GADest(GlobalArray<ndim,T>& ga, coords_range<ndim> r);
		public:
			~GADest();
			GADest(const GADest<ndim,T>& gad);
			GADest& operator=(const GADest<ndim,T>& gad) = delete;
			template<typename S>
			GADest<ndim,T>& operator=(const S& src);
		};

		template<int ndim, typename T>
		class GARef {
			const GlobalArray<ndim,T>& ga;
		public:
			GARef(const GlobalArray<ndim,T>& ga);
			~GARef();
			INLINE operator const GlobalArray<ndim,T>&() const;
			INLINE const Domain<ndim>& domain() const;
			INLINE coords_range<ndim> region() const;
		};

		// Inline function implementations

		template<int ndim, typename T>
		inline const Domain<ndim>& GlobalArray<ndim,T>::domain() const {
			return *dom;
		}

		template<int ndim, typename T>
		inline coords<ndim> GlobalArray<ndim,T>::ghost_width() const {
			return gw;
		}

		template<int ndim, typename T>
		inline bool GlobalArray<ndim,T>::ghost_corners() const {
			return gc;
		}

		template<int ndim, typename T>
		inline GlobalArray<ndim,T>::operator bool() const {
			return dom != 0;
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::region() const {
			return domain().total();
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::inner() const {
			auto r = domain().local();
                        r.lower += ghost_width();
                        r.upper -= ghost_width();
                        return r;
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::outer_front(int d) const {
                        assert(d >= 0 && d < ndim);
			auto r = domain().local();
                        r.upper[d] = r.lower[d] + ghost_width()[d];
                        return r;
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::outer_back(int d) const {
                        assert(d >= 0 && d < ndim);
			auto r = domain().local();
                        r.lower[d] = r.upper[d] - ghost_width()[d];
                        return r;
		}

		template<int ndim, typename T>
		inline T& GlobalArray<ndim,T>::da(coords<ndim> i) const {
			return ptr[(i + gw).offset(ld)];
		}

		template<int ndim, typename T>
		inline GARef<ndim,T>::operator const GlobalArray<ndim,T>&() const {
			return ga;
		}

		template<int ndim, typename T>
		inline const Domain<ndim>& GARef<ndim,T>::domain() const {
			return ga.domain();
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GARef<ndim,T>::region() const {
			return ga.region();
		}

		// Generic members

		template<int ndim, typename T> template<typename S>
		GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator=(const S& src) {
			region(domain().total()) = src;
			return *this;
		}

		template<int ndim, typename T> template<typename S>
		GADest<ndim,T>& GADest<ndim,T>::operator=(const S& src) {
			static_assert(S::number_of_dimensions == ndim, "source dimensionality");
			assert(ga.domain() == src.domain());
			assert(src.region().contains(r));
			Access<ndim,T> d(ga);
			const typename S::accessor s(src);
			ga.domain().for_each(r, [&d, &s](coords<ndim> i){
				d(i) = s(i);
			});
			return *this;
		}
	}

}

#endif
