#ifndef __SHARK_REGION_HPP
#define __SHARK_REGION_HPP

#include "common.hpp"
#include "domain.hpp"
#include "coords_range.hpp"

namespace shark {

	namespace ndim {
		template<int ndim>
		class Region {
		private:
			const Domain<ndim>& domain;
			const coords_range<ndim> r;

		public:
			Region(const Domain<ndim>& domain);
			Region(const Domain<ndim>& domain, coords_range<ndim> r);
			~Region();

			bool hasData() const;

			template<typename R, typename Func>
			R internal_add_reduce(const R& zero, const Func& f) const;
		
			/**
			* Apply elemental function onto the elements of range
			*/
			template<typename Func>
			void for_each(const Func& f) const;

		};

		// Inline Region member functions

		template<int ndim> template<typename Func>
		void Region<ndim>::for_each(const Func& f) const {
#if defined(SHARK_SER_SCHED)
			r.for_each(f);
#else
#error "No scheduler for_each"
#endif
		}

		template<int ndim> template<typename R, typename Func>
		R Region<ndim>::internal_add_reduce(const R& identity, const Func& f) const {
#if defined(SHARK_SER_SCHED)
			R local_sum = identity;
			r.add(f, local_sum);
			return local_sum;
#else
#error "No scheduler internal_add_reduce"
#endif
		}
		

	}

}

#endif
