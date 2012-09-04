#ifndef __SHARK_COORDS_HPP
#define __SHARK_COORDS_HPP

#include <array>                       // std::array
#include <cstddef>                     // std::size_t
#include <ostream>                     // std::ostream
#include "common.hpp"
#include "vec.hpp"

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
			INLINE coords<ndim> operator-() const;
			INLINE coords<ndim> operator+(const coords<ndim>& other) const;
			INLINE coords<ndim> operator-(const coords<ndim>& other) const;

			INLINE coords<ndim>& operator<<=(unsigned short w);
			INLINE coords<ndim>& operator>>=(unsigned short w);
			INLINE coords<ndim> operator<<(unsigned short w) const;
			INLINE coords<ndim> operator>>(unsigned short w) const;

			INLINE coord offset(const coords<ndim+1>& ld) const;
			INLINE coord offset(const std::array<std::size_t,ndim-1>& ld) const;

			INLINE vec<ndim,coord> to_vec() const;

			INLINE static coords<ndim> one();
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
			return seq<0,ndim>::all_of([this,&other](int d) { return val[d] == other.val[d]; });
		}

		template<int ndim>
		inline bool coords<ndim>::operator!=(const coords<ndim>& other) const {
			return !(*this == other);
		}

		template<int ndim>
		inline bool coords<ndim>::operator<(const coords<ndim>& other) const {
			return seq<0,ndim>::all_of([this,&other](int d) { return val[d] < other.val[d]; });
		}

		template<int ndim>
		inline bool coords<ndim>::operator<=(const coords<ndim>& other) const {
			return seq<0,ndim>::all_of([this,&other](int d) { return val[d] <= other.val[d]; });
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
			seq<0,ndim>::for_each([this,&other](int d) { val[d] += other.val[d]; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator-() const {
			coords<ndim> r;
			seq<0,ndim>::for_each([this,&r](int d) { r.val[d] = -val[d]; });
			return r;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator+(const coords<ndim>& other) const {
			return coords<ndim>(*this) += other;
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator-=(const coords<ndim>& other) {
			seq<0,ndim>::for_each([this,&other](int d) { val[d] -= other.val[d]; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::operator-(const coords<ndim>& other) const {
			return coords<ndim>(*this) -= other;
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator<<=(unsigned short w) {
			seq<0,ndim>::for_each([this,w](int d) { val[d] <<= w; });
			return *this;
		}

		template<int ndim>
		inline coords<ndim>& coords<ndim>::operator>>=(unsigned short w) {
			seq<0,ndim>::for_each([this,w](int d) { val[d] >>= w; });
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
			return seq<0,ndim>::sum([this,&ld](int d) { return val[d] * ld[d+1]; }, coord(0));
		}

		template<int ndim>
		inline coord coords<ndim>::offset(const std::array<std::size_t,ndim-1>& ld) const {
			return seq<0,ndim-1>::sum([this,&ld](int d) { return val[d] * static_cast<coord>(ld[d]); }, val[ndim-1]);
		}

		template<int ndim>
		inline vec<ndim,coord> coords<ndim>::to_vec() const {
			vec<ndim,coord> r;
			seq<0,ndim>::for_each([this,&r](int d) { r[d] = (*this)[d]; });
			return r;
		}

		template<int ndim>
		inline coords<ndim> coords<ndim>::one() {
			coords<ndim> r;
			seq<0,ndim>::for_each([&r](int d) { r[d] = coord(1); });
			return r;
		}

	}

}

#endif
