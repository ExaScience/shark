#ifndef __SHARK_ACCESS_HPP
#define __SHARK_ACCESS_HPP

#include "common.hpp"
#include "coords.hpp"
#include "domain.hpp"
#include "globalarray.hpp"

namespace shark {

	namespace ndim {

		/**
		 * Direct access to the local region of a GlobalArray.
		 * Local access cannot overlap with RMA operations.
		 */
		template<int ndim, typename T>
		class Access {
			friend class Region<ndim>;

		private:
			const GlobalArray<ndim,T>* ga;
			void release();
			void reset();

		public:
			/**
			 * Construct active access to the local region of a GlobalArray (local).
			 * @param ga the global array
			 */
			Access(const GlobalArray<ndim,T>& ga);
			Access(GlobalArray<ndim,T>& ga);

			/**
			 * Returns whether this is an active access.
			 */
			INLINE /*explicit*/ operator bool() const;

			/**
			 * Create an inactive access.
			 */
			Access();

			// No copying semantics
			Access(const Access<ndim,T>&) = delete;
			Access& operator=(const Access<ndim,T>&) = delete;

			// Move semantics
			Access(Access<ndim,T>&& other);
			Access& operator=(Access<ndim,T>&& other);

			/**
			 * Terminate an Access.
			 * If active, this will release the access on the GlobalArray.
			 */
			~Access();

			/**
			 * Access an element
			 * @param i the n-dimensional index of the element to access
			 */
			INLINE const T& operator()(coords<ndim> i) const;
			INLINE T& operator()(coords<ndim> i);

			/**
			 * Access an element using a local index
			 * @param i the n-dimensional local index of the element to access
			 */
			INLINE const T& getLocal(coords<ndim> i) const;
			INLINE T& getLocal(coords<ndim> i);
		};

		// Inline Access member functions
		
		template<int ndim, typename T>
		inline Access<ndim,T>::operator bool() const {
			return ga != 0;
		}

		template<int ndim,typename T>
		inline T& Access<ndim,T>::operator()(coords<ndim> i) {
			return ga->access(i - ga->domain().local().lower);
		}

		template<int ndim,typename T>
		inline const T& Access<ndim,T>::operator()(coords<ndim> i) const {
			return ga->access(i - ga->domain().local().lower);
		}

		template<int ndim,typename T>
		inline T& Access<ndim,T>::getLocal(coords<ndim> i) {
			return ga->access(i);
		}

		template<int ndim,typename T>
		inline const T& Access<ndim,T>::getLocal(coords<ndim> i) const {
			return ga->access(i);
		}

	}

}

#endif
