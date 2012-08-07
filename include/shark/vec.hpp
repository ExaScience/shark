#ifndef __SHARK_VEC_HPP
#define __SHARK_VEC_HPP

#include "common.hpp"

namespace shark {

	namespace ndim {

		template<int ndim,typename T>
		struct vec {
			T val[ndim];

			// No explicit construct/copy/destroy for aggregate type.

			INLINE T& operator[](int d);
			INLINE const T& operator[](int d) const;

			INLINE bool operator==(const vec<ndim,T>& other) const;
			INLINE bool operator!=(const vec<ndim,T>& other) const;

			INLINE vec<ndim,T>& operator+=(const vec<ndim,T>& other);
			INLINE vec<ndim,T> operator+(const vec<ndim,T>& other) const;
		};

		template<int ndim,typename T>
		inline T& vec<ndim,T>::operator[](int d) {
			return val[d];
		}

		template<int ndim,typename T>
		inline const T& vec<ndim,T>::operator[](int d) const {
			return val[d];
		}

		template<int ndim,typename T>
		inline bool vec<ndim,T>::operator==(const vec<ndim,T>& other) const {
			return seq<0,ndim>::all_of([this,&other](int d) { return this->val[d] == other.val[d]; } );
		}
		
		template<int ndim,typename T>
		inline bool vec<ndim,T>::operator!=(const vec<ndim,T>& other) const {
			return !(*this == other);
		}

		template<int ndim,typename T>
		inline vec<ndim,T>& vec<ndim,T>::operator+=(const vec<ndim,T>& other) {
			seq<0,ndim>::for_each([this,&other](int d) { this->val[d] += other.val[d]; });
			return *this;
		}

		template<int ndim,typename T>
		inline vec<ndim,T> vec<ndim,T>::operator+(const vec<ndim,T>& other) const {
			return vec<ndim,T>(*this) += other;
		}

	}

}
#endif
