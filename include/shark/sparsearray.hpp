/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_SPARSEARRAY_HPP
#define __SHARK_SPARSEARRAY_HPP

#include <algorithm>                   // std::min
#include <array>                       // std::array
#include <cassert>                     // assert
#include <cstddef>                     // std::size_t
#include <vector>                      // std::vector
#include "common.hpp"
#include "coords.hpp"
#include "coords_range.hpp"

namespace shark {

	namespace ndim {

		/**
		 * A local data structure to store values for specific positions in a global domain.
		 * The data structure is sparse: only selection of the positions is available.
		 * This is an auxiliary data structure for gather/scatter operations.
		 */
		template<int ndim,typename T>
		class SparseArray {
			template<int,typename> friend class SparseArray;
			friend class GlobalArray<ndim,T>;

		private:
			const Domain<ndim>& dom;
			const unsigned short cw;
			const coords<ndim+1> ld;
			const coords<ndim+1> cld;
			T* const ptr;
			std::vector<bool> coarse;

			static unsigned short full_cw(const Domain<ndim>& dom);
			static coords<ndim+1> stride(const Domain<ndim>& dom);
			static coords<ndim+1> coarse_stride(const Domain<ndim>& dom, unsigned short cw);

			void init();

			template<typename Func>
			void iter(const Func& f) const;

#ifdef SHARK_MPI_COMM
			void proc_ranges(std::vector<coords_range<ndim>> local[], std::vector<coords_range<ndim>> global[]) const;
#endif

		public:
			/**
			 * The domain of the sparse array
			 */
			INLINE const Domain<ndim>& domain() const;

			/**
			 * Create a sparse array for a given domain
			 * @param dom the domain
			 * @param cw the coarsening width (the number of coordinate digits to drop) (default: full coarsen)
			 */
			SparseArray(const Domain<ndim>& dom, unsigned short cw);
			SparseArray(const Domain<ndim>& dom);

			/**
			 * Copy a sparse array
			 * @param other the sparse array to copy from
			 */
			SparseArray(const SparseArray<ndim,T>& other);
			SparseArray<ndim,T>& operator=(const SparseArray<ndim,T>& other);

			/**
			 * Create a sparse array from another sparse array of unrelated type
			 * The new sparse array will have (at least) the same positions available as the original one
			 * (but the corresponding values will be undefined)
			 * @param other the sparse array to copy from
			 */
			template<typename U>
			SparseArray(const SparseArray<ndim,U>& other);

			/**
			 * Destruct a sparse array
			 */
			~SparseArray();

			/**
			 * Empty this sparse array; no positions will be available
			 */
			void clear();

			/**
			 * Store a value in a given position. This also marks the position as available. If this
			 * sparse array is used as destination of a gather operation, the position will be filled
			 * with a remote value.
			 * @param i the coordinates of the position
			 * @param v the value
			 */
			INLINE void set(coords<ndim> i, T v = T());

			/**
			 * Increase a value in a given position. This will mark the position as available.
			 * @param i the coordinates of the position
			 * @param v the value to increase by
			 * @see set(coords<ndim>, T)
			 */
			INLINE void increment(coords<ndim> i, T v);

			/**
			 * Retrieve the value in a given position. This assumes that the position is available, otherwise
			 * the result is undefined.
			 * @param i the coordinates of the position
			 */
			INLINE T get(coords<ndim> i) const;

			/**
			 * Accumulate values of another SparseArray.
			 */
			void operator+=(const SparseArray<ndim,T>& other);
		};

		// Inline SparseArray members
		template<int ndim,typename T>
		inline const Domain<ndim>& SparseArray<ndim,T>::domain() const {
			return dom;
		}

		template<int ndim,typename T>
		inline void SparseArray<ndim,T>::set(coords<ndim> i, T v) {
			coarse[(i >> cw).offset(cld)] = true;
			ptr[i.offset(ld)] = v;
		}

		template<int ndim,typename T>
		inline void SparseArray<ndim,T>::increment(coords<ndim> i, T v) {
			coarse[(i >> cw).offset(cld)] = true;
			ptr[i.offset(ld)] += v;
		}

		template<int ndim,typename T>
		inline T SparseArray<ndim,T>::get(coords<ndim> i) const {
			assert(coarse[(i >> cw).offset(cld)]);
			return ptr[i.offset(ld)];
		}

		template<int ndim,typename T> template<typename U>
		SparseArray<ndim,T>::SparseArray(const SparseArray<ndim,U>& other):
			dom(other.dom), cw(other.cw), ld(other.ld), cld(other.cld), ptr(new T[ld[0]]), coarse(other.coarse) {
			init();
		}

		template<int ndim,typename T> template<typename Func>
		void SparseArray<ndim,T>::iter(const Func& f) const {
			for(coord k = 0; k < cld[0]; k++)
				if(coarse[k]) {
					coords_range<ndim> r;
					coord kk = k;
					for(int d = 0; d < ndim; d++) {
						coord fk = kk / cld[d+1];
						r.lower[d] = fk << cw;
						r.upper[d] = std::min((fk+1) << cw, dom.n[d]);
						kk = kk % cld[d+1];
					}
					f(r);
				}
		}

	}

}

#endif
