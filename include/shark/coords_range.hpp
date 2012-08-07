
#ifndef __SHARK_COORDS_RANGE_HPP
#define __SHARK_COORDS_RANGE_HPP

#include <ostream>                     // std::ostream
#include "common.hpp"
#include "coords.hpp"

namespace shark {
	
	namespace ndim {

		/**
		 * An n-dimensional coordinate range.
		 */
		template<int ndim>
		struct coords_range {
			coords<ndim> lower;
			coords<ndim> upper;

			// No explicit construct/copy/destroy for aggregate type.
			
		private:
			template<int d, typename Func>
			INLINE typename std::enable_if<d < ndim>::type for_eachd(const Func& f, coords<ndim>& i) const;
			template<int d, typename Func>
			INLINE typename std::enable_if<d == ndim>::type for_eachd(const Func& f, coords<ndim>& i) const;

		public:
			template<typename Func>
			INLINE void for_each(const Func& f) const;
			template<typename Func, typename Acc>
			INLINE void add(const Func& f, Acc& a) const;

			coords_range<ndim> overlap(coords_range<ndim> other) const;
			bool contains(coords<ndim> i) const;
			bool contains(coords_range<ndim> other) const;
			INLINE coord size() const;
		};

		template<int ndim>
		std::ostream& operator<<(std::ostream& out, const coords_range<ndim>& r);

		// Inline function implementations

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d < ndim>::type coords_range<ndim>::for_eachd(const Func& f, coords<ndim>& i) const {
			for(i[d] = lower[d]; i[d] < upper[d]; i[d]++)
				for_eachd<d+1>(f, i);
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d == ndim>::type coords_range<ndim>::for_eachd(const Func& f, coords<ndim>& i) const {
			f(i);
		}

		template<int ndim> template<typename Func>
		inline void coords_range<ndim>::for_each(const Func& f) const {
			coords<ndim> i;
			for_eachd<0>(f, i);
		}

		template<int ndim> template<typename Func, typename Acc>
		inline void coords_range<ndim>::add(const Func& f, Acc& a) const {
			for_each([&f, &a](coords<ndim> i) { a += f(i); });
		}

		template<int ndim>
		inline coord coords_range<ndim>::size() const {
			return seq<0,ndim>::product([this](int d) { return upper[d] - lower[d]; }, coord(1));
		}

	}

}

#endif
