#ifndef __SHARK_VEC_HPP
#define __SHARK_VEC_HPP

#include <utility>        // std::declval
#include "abs.hpp"
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

			template<typename S>
			INLINE vec<ndim,T>& operator+=(const vec<ndim,S>& other);
			template<typename S>
			INLINE vec<ndim,T>& operator-=(const vec<ndim,S>& other);
			template<typename S>
			INLINE vec<ndim,T>& operator*=(const vec<ndim,S>& other);
			template<typename S>
			INLINE vec<ndim,T>& operator/=(const vec<ndim,S>& other);

			INLINE T min() const;
			INLINE T max() const;

			INLINE static vec<ndim,T> one();
		};

		template<int ndim, typename T>
		INLINE vec<ndim,decltype(-std::declval<T>())> operator-(const vec<ndim,T>& v);

		template<int ndim, typename T, typename S>
		INLINE vec<ndim,decltype(std::declval<T>() + std::declval<S>())> operator+(const vec<ndim,T>& v1, const vec<ndim,S>& v2);

		template<int ndim, typename T, typename S>
		INLINE vec<ndim,decltype(std::declval<T>() - std::declval<S>())> operator-(const vec<ndim,T>& v1, const vec<ndim,S>& v2);

		// Declaring operator* triggers ambiguous overload
		// http://llvm.org/bugs/show_bug.cgi?id=13730
		//template<int ndim, typename T, typename S>
		//INLINE vec<ndim,decltype(std::declval<T>() * std::declval<S>())> operator*(const vec<ndim,T>& v1, const vec<ndim,S>& v2);

		//template<int ndim, typename T, typename S>
		//INLINE vec<ndim,decltype(std::declval<T>() * std::declval<S>())> operator*(const vec<ndim,T>& v, const S& a);

		//template<int ndim, typename T, typename S>
		//INLINE vec<ndim,decltype(std::declval<T>() * std::declval<S>())> operator*(const T& a, const vec<ndim,S>& v);

		template<int ndim, typename T, typename S>
		INLINE vec<ndim,decltype(std::declval<T>() / std::declval<S>())> operator/(const vec<ndim,T>& v1, const vec<ndim,S>& v2);

		//Declaring abs triggers ambiguous overlad 
		//http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43054
		//template<int ndim, typename T>
		//INLINE vec<ndim,decltype(abs(std::declval<T>()))> abs(const vec<ndim,T>& v);

		// Inline / generic members

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
			return seq<0,ndim>::all_of([this,&other](int d) { return (*this)[d] == other[d]; } );
		}
		
		template<int ndim,typename T>
		inline bool vec<ndim,T>::operator!=(const vec<ndim,T>& other) const {
			return !(*this == other);
		}

		template<int ndim,typename T> template<typename S>
		inline vec<ndim,T>& vec<ndim,T>::operator+=(const vec<ndim,S>& other) {
			seq<0,ndim>::for_each([this,&other](int d) { (*this)[d] += other[d]; });
			return *this;
		}

		template<int ndim,typename T> template<typename S>
		inline vec<ndim,T>& vec<ndim,T>::operator-=(const vec<ndim,S>& other) {
			seq<0,ndim>::for_each([this,&other](int d) { (*this)[d] -= other[d]; });
			return *this;
		}
		template<int ndim,typename T> template<typename S>
		inline vec<ndim,T>& vec<ndim,T>::operator*=(const vec<ndim,S>& other) {
			seq<0,ndim>::for_each([this,&other](int d) { (*this)[d] *= other[d]; });
			return *this;
		}
		template<int ndim,typename T> template<typename S>
		inline vec<ndim,T>& vec<ndim,T>::operator/=(const vec<ndim,S>& other) {
			seq<0,ndim>::for_each([this,&other](int d) { (*this)[d] /= other[d]; });
			return *this;
		}

		template<int ndim,typename T>
		inline T vec<ndim,T>::min() const {
			T cand((*this)[0]);
			seq<1,ndim-1>::for_each([this,&cand](int d) {
				if((*this)[d] < cand)
					cand = (*this)[d];
			});
			return cand;
		}

		template<int ndim,typename T>
		inline T vec<ndim,T>::max() const {
			T cand((*this)[0]);
			seq<1,ndim-1>::for_each([this,&cand](int d) {
				if((*this)[d] > cand)
					cand = (*this)[d];
			});
			return cand;
		}

		template<int ndim, typename T>
		inline vec<ndim,T> vec<ndim,T>::one() {
			vec<ndim,T> r;
			seq<0,ndim>::for_each([&r](int d) { r[d] = T(1); });
			return r;
		}

		template<int ndim, typename T>
		inline vec<ndim,decltype(-std::declval<T>())> operator-(const vec<ndim,T>& v) {
			vec<ndim,decltype(-std::declval<T>())> r;
			seq<0,ndim>::for_each([&r,&v](int d) { r[d] = -v[d]; });
			return r;
		}

		template<int ndim, typename T, typename S>
		inline vec<ndim,decltype(std::declval<T>() + std::declval<S>())> operator+(const vec<ndim,T>& v1, const vec<ndim,S>& v2) {
			vec<ndim,decltype(std::declval<T>() + std::declval<S>())> r;
			seq<0,ndim>::for_each([&r,&v1,&v2](int d) { r[d] = v1[d] + v2[d]; });
			return r;
		}

		template<int ndim, typename T, typename S>
		inline vec<ndim,decltype(std::declval<T>() - std::declval<S>())> operator-(const vec<ndim,T>& v1, const vec<ndim,S>& v2) {
			vec<ndim,decltype(std::declval<T>() - std::declval<S>())> r;
			seq<0,ndim>::for_each([&r,&v1,&v2](int d) { r[d] = v1[d] - v2[d]; });
			return r;
		}

		template<int ndim, typename T, typename S>
		inline vec<ndim,decltype(std::declval<T>() * std::declval<S>())> operator*(const vec<ndim,T>& v1, const vec<ndim,S>& v2) {
			vec<ndim,decltype(std::declval<T>() * std::declval<S>())> r;
			seq<0,ndim>::for_each([&r,&v1,&v2](int d) { r[d] = v1[d] * v2[d]; });
			return r;
		}

		template<int ndim, typename T, typename S>
		inline vec<ndim,decltype(std::declval<T>() * std::declval<S>())> operator*(const vec<ndim,T>& v, const S& a) {
			vec<ndim,decltype(std::declval<T>() * std::declval<S>())> r;
			seq<0,ndim>::for_each([&r,&v,&a](int d) { r[d] = v[d] * a; });
			return r;
		}

		template<int ndim, typename T, typename S>
		inline vec<ndim,decltype(std::declval<T>() * std::declval<S>())> operator*(const T& a, const vec<ndim,S>& v) {
			vec<ndim,decltype(std::declval<T>() * std::declval<S>())> r;
			seq<0,ndim>::for_each([&r,&a,&v](int d) { r[d] = a * v[d]; });
			return r;
		}

		template<int ndim, typename T, typename S>
		inline vec<ndim,decltype(std::declval<T>() / std::declval<S>())> operator/(const vec<ndim,T>& v1, const vec<ndim,S>& v2) {
			vec<ndim,decltype(std::declval<T>() / std::declval<S>())> r;
			seq<0,ndim>::for_each([&r,&v1,&v2](int d) { r[d] = v1[d] / v2[d]; });
			return r;
		}

		template<int ndim, typename T>
		inline vec<ndim,decltype(abs(std::declval<T>()))> abs(const vec<ndim,T>& v) {
			vec<ndim,decltype(abs(std::declval<T>()))> r;
			seq<0,ndim>::for_each([&r,&v](int d) { r[d] = abs(v[d]); });
			return r;
		}

	}

}
#endif
