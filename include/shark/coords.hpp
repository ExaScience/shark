#ifndef __SHARK_COORDS_HPP
#define __SHARK_COORDS_HPP

#include <ostream>                     // std::ostream
#include "common.hpp"

namespace shark {

	namespace ndim {

		/**
		 * An n-dimensional coordinate pair.
		 */
		template<int ndim>
		struct coords {
			coord val[ndim];

			// No explicit construct/copy/destroy for aggregate type.

			INLINE coord& operator[](int d);
			INLINE const coord& operator[](int d) const;

			INLINE bool operator==(const coords<ndim>& other) const;
			INLINE bool operator!=(const coords<ndim>& other) const;
			INLINE bool operator< (const coords<ndim>& other) const;
			INLINE bool operator<=(const coords<ndim>& other) const;
			INLINE bool operator> (const coords<ndim>& other) const;
			INLINE bool operator>=(const coords<ndim>& other) const;

			INLINE coords<ndim>& operator+=(const coords<ndim>& other);
			INLINE coords<ndim>& operator-=(const coords<ndim>& other);
			INLINE coords<ndim> operator+(const coords<ndim>& other) const;
			INLINE coords<ndim> operator-(const coords<ndim>& other) const;

			INLINE coords<ndim>& operator<<=(unsigned short w);
			INLINE coords<ndim>& operator>>=(unsigned short w);
			INLINE coords<ndim> operator<<(unsigned short w) const;
			INLINE coords<ndim> operator>>(unsigned short w) const;

			INLINE coord offset(const coords<ndim+1>& ld) const;
		};

		template<int ndim>
		std::ostream& operator<<(std::ostream& out, const coords<ndim>& i);


		// Inline function implementations

		template<int ndim>
		inline coord& coords<ndim>::operator[](int d) {
			return val[d];
		}

		template<int ndim>
		inline const coord& coords<ndim>::operator[](int d) const {
			return val[d];
		}

		template<int ndim>
		inline bool coords<ndim>::operator==(const coords<ndim>& other) const {
			return all_of<0,ndim>([this,&other](int d) { return val[d] == other.val[d]; });
		}

		template<int ndim>
		inline bool coords<ndim>::operator!=(const coords<ndim>& other) const {
			return !(*this == other);
		}

		template<int ndim>
		inline bool coords<ndim>::operator<(const coords<ndim>& other) const {
			return all_of<0,ndim>([this,&other](int d) { return val[d] < other.val[d]; });
		}

		template<int ndim>
		inline bool coords<ndim>::operator<=(const coords<ndim>& other) const {
			return all_of<0,ndim>([this,&other](int d) { return val[d] <= other.val[d]; });
		}

		template<int ndim>
		inline bool coords<ndim>::operator>(const coords<ndim>& other) const {
			return !(*this <= other);
		}

		template<int ndim>
		inline bool coords<ndim>::operator>=(const coords<ndim>& other) const {
			return !(*this < other);
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator+=(const coords<ndim>& other) {
			for_each<0,ndim>([this,&other](int d) { val[d] += other.val[d]; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator+(const coords<ndim>& other) const {
			return coords<ndim>(*this) += other;
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator-=(const coords<ndim>& other) {
			for_each<0,ndim>([this,&other](int d) { val[d] -= other.val[d]; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator-(const coords<ndim>& other) const {
			return coords<ndim>(*this) -= other;
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator<<=(unsigned short w) {
			for_each<0,ndim>([this,w](int d) { val[d] <<= w; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator>>=(unsigned short w) {
			for_each<0,ndim>([this,w](int d) { val[d] >>= w; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator<<(unsigned short w) const {
			return coords<ndim>(*this) <<= w;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator>>(unsigned short w) const {
			return coords<ndim>(*this) >>= w;
		}

		template<int ndim>
		inline coord coords<ndim>::offset(const coords<ndim+1>& ld) const {
			return sum<0,ndim>([this,&ld](int d) { return val[d] * ld[d+1]; }, coord(0));
		}

	}

}

#endif
