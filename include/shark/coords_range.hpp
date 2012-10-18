
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

			coords_range<ndim> overlap(coords_range<ndim> other) const;
			bool contains(coords<ndim> i) const;
			bool contains(coords_range<ndim> other) const;
			INLINE coord size() const;
			INLINE coords<ndim> periodic_equiv(coords<ndim> i) const;
			INLINE coords<ndim> adj(int d, int shift) const;
			coords<ndim+1> stride(coords<ndim> bw = coords<ndim>()) const;
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

		template<int ndim>
		inline coord coords_range<ndim>::size() const {
			return seq<0,ndim>::product([this](int d) { return upper[d] - lower[d]; }, coord(1));
		}

		template<int ndim>
		inline coords<ndim> coords_range<ndim>::periodic_equiv(coords<ndim> ii) const {
			seq<0,ndim>::for_each([this,&ii](int d) {
				coord nd = upper[d] - lower[d];
				if(ii[d] < lower[d]) {
					do {
						ii[d] += nd;
					} while(ii[d] < lower[d]);
				} else if(ii[d] >= upper[d]) {
					do {
						ii[d] -= nd;
					} while(ii[d] >= upper[d]);
				}
			});
			return ii;
		}

		template<int ndim>
		inline coords<ndim> coords_range<ndim>::adj(int d, int shift) const {
			coords<ndim> adj = lower;
			adj[d] = lower[d] + shift * (upper[d] - lower[d]);
			return adj;
		}

	}

}

#endif
