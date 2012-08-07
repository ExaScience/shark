#ifndef __SHARK_GLOBALARRAY_HPP
#define __SHARK_GLOBALARRAY_HPP

#include <array>                       // std::array
#include <memory>                      // std::unique_ptr
#include "common.hpp"
#include "coords.hpp"
#include "coords_range.hpp"

namespace shark {

	namespace ndim {

		template<int>
		class GAImpl;

		/**
		 * A global array that is partitioned across a process group according to a domain
		 */
		template<int ndim,typename T>
		class GlobalArray {
			friend class Access<ndim,T>;
		
		public:
			typedef Access<ndim,T> accessor;
			static const /*constexpr*/ int number_of_dimensions = ndim;

		private:
			const Domain<ndim>* dom;
			coords<ndim> gw;
			bool gc;
			T* ptr;

			std::unique_ptr<GAImpl<ndim>> impl;
			coords<ndim+1> ld;

			mutable int lc;

			static void allocate(const Domain<ndim>& domain, coords<ndim> ghost_width, bool ghost_corners, T **ptr, GAImpl<ndim>& impl, coords<ndim+1>& ld);
			void deallocate();
			void reset();

			static coords<ndim+1> stride(coords<ndim> count, coords<ndim> ghost_width = coords<ndim>());

			INLINE T& access(coords<ndim> i) const;

			class RMAOp;

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
			GlobalArray(const Domain<ndim>& domain, coords<ndim> ghost_width = coords<ndim>(), bool ghost_corners = false);

			/**
			 * Construct a GlobalArray as a copy of another (collective).
			 * The new array will have the same domain and ghost configuration as the old one.
			 * @param other the array to copy from, \verbatim hasSameShape(other) \endverbatim
			 * @param copy whether to copy the data from the other array
			 */
			GlobalArray(const GlobalArray<ndim,T>& other, bool copy);

			/**
			 * Destruct a GlobalArray (collective).
			 * If active, the memory of the global array is released
			 */
			~GlobalArray();

			// No copy semantics
			GlobalArray(const GlobalArray<ndim,T>&) = delete;
			GlobalArray<ndim,T>& operator=(const GlobalArray<ndim,T>&) = delete;

			// Move semantics
			GlobalArray(GlobalArray<ndim,T>&& other);
			GlobalArray<ndim,T>& operator=(GlobalArray<ndim,T>&& other);

			/**
			 * Redistribute the GlobalArray to have a Domain.
			 * @param domain the new domain, \verbatim domain().equiv(domain) \enverbatim
			 */
			void reshape(const Domain<ndim>& domain);

			/**
			 * Check whether some other GlobalArray has the same domain as this one (local).
			 * @param other the global array to check against
			 */
			bool hasSameShape(const GlobalArray<ndim,T>& other) const;

			/**
			 * Update the ghost cells with data from their original locations (collective).
			 * RMA operations cannot overlap with local access.
			 */
			void update() const;

			/**
			 * Get remote range (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to retrieve
			 * @param ld the strides to use for buf
			 * @param buf the target buffer
			 */
			void get(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, T* buf) const;

			/**
			 * Put remote range (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to retrieve
			 * @param ld the strides to use for buf
			 * @param buf the target buffer
			 */
			void put(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, const T* buf);

			/**
			 * Increase remote range with local values (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to update
			 * @param ld the strides to use for buf
			 * @param buf the source buffer
			 */
			void accumulate(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, const T* buf);

			/**
			 * Selective gather of remote values (one-sided). This fills all available positions in a
			 * local sparse array with their values from this global array.
			 * RMA operations cannot overlap with local access.
			 * @param sa the sparse array to fill
			 */
			void gather(SparseArray<ndim,T>& sa) const;

			/**
			 * Selective scatter to remote values (one-sided). This sends out all available values of a
			 * local sparse array to the respective positions in this global array. The values are accumulated
			 * with the original value at the destination.
			 * RMA operations cannot overlap with local access.
			 * @param sa the sparse array whose values to send out
			 */
			void scatterAcc(const SparseArray<ndim,T>& sa);
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
			return ptr != 0;
		}

		template<int ndim, typename T>
		inline T& GlobalArray<ndim,T>::access(coords<ndim> i) const {
			return ptr[(i + gw).offset(ld)];
		}

	}

}

#endif
