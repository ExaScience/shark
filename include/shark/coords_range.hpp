
#ifndef __SHARK_COORDS_RANGE_HPP
#define __SHARK_COORDS_RANGE_HPP

#include <ostream>                     // std::ostream
#include "common.hpp"
#include "coords.hpp"

#ifdef __INTEL_COMPILER
#define SHARK_VECTOR_PRAGMA
#endif

namespace shark {
	
	namespace ndim {

		/**
		 * An n-dimensional coordinate range.
		 */
		template<int ndim>
		struct coords_range {
			friend class Domain<ndim>;

			coords<ndim> lower;
			coords<ndim> upper;

			// No explicit construct/copy/destroy for aggregate type.
			
		private:
			template<int d, typename Func>
			INLINE typename std::enable_if<d < ndim-1>::type for_eachd(const Func& f, coords<ndim>& i) const;
			template<int d, typename Func>
			INLINE typename std::enable_if<d == ndim-1>::type for_eachd(const Func& f, coords<ndim>& i) const;

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
		inline typename std::enable_if<d < ndim-1>::type coords_range<ndim>::for_eachd(const Func& f, coords<ndim>& i) const {
			for(i[d] = lower[d]; i[d] < upper[d]; i[d]++)
				for_eachd<d+1>(f, i);
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d == ndim-1>::type coords_range<ndim>::for_eachd(const Func& f, coords<ndim>& i) const {
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
			for(coord j = lower[d]; j < upper[d]; j++) {
				i[d] = j;
				f(i);
			}
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

#if defined(SHARK_OMP_SCHED)
		template<int>
		struct omp_coords_range;

		template<>
		struct omp_coords_range<1> {
			coords_range<1> r;

			template<typename Func>
			INLINE void for_each(const Func& f) const;

			template<typename T, typename Func>
			INLINE typename std::enable_if<std::is_scalar<T>::value,T>::type
			internal_sum(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			INLINE typename std::enable_if<!std::is_scalar<T>::value,T>::type
			internal_sum(const T& zero, const Func& f) const;
		};

		template<typename Func>
		inline void omp_coords_range<1>::for_each(const Func& f) const {
#pragma omp parallel for schedule(static)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
			for(coord j = r.lower[0]; j < r.upper[0]; j++) {
				coords<1> i;
				i[0] = j;
				f(i);
			}
		}

		template<typename T, typename Func>
		inline typename std::enable_if<std::is_scalar<T>::value,T>::type
		omp_coords_range<1>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
#pragma omp parallel for schedule(static) reduction(+: sum)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
			for(coord j = r.lower[0]; j < r.upper[0]; j++) {
				coords<1> i;
				i[0] = j;
				f(sum, i);
			}
			return sum;
		}

		template<typename T, typename Func>
		inline typename std::enable_if<!std::is_scalar<T>::value,T>::type
		omp_coords_range<1>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
#pragma omp parallel shared(sum)
			{
				T local_sum(zero);
#pragma omp for schedule(static)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
				for(coord j = r.lower[0]; j < r.upper[0]; j++) {
					coords<1> i;
					i[0] = j;
					f(local_sum, i);
				}
#pragma omp critical
				sum += local_sum;
			}
			return sum;
		}

		template<>
		struct omp_coords_range<2> {
			coords_range<2> r;

			template<typename Func>
			INLINE void for_each(const Func& f) const;

			template<typename T, typename Func>
			INLINE typename std::enable_if<std::is_scalar<T>::value,T>::type
			internal_sum(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			INLINE typename std::enable_if<!std::is_scalar<T>::value,T>::type
			internal_sum(const T& zero, const Func& f) const;
		};

		template<typename Func>
		inline void omp_coords_range<2>::for_each(const Func& f) const {
#pragma omp parallel for schedule(static) collapse(2)
			for(coord j = r.lower[0]; j < r.upper[0]; j++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
				for(coord k = r.lower[1]; k < r.upper[1]; k++) {
					coords<2> i;
					i[0] = j;
					i[1] = k;
					f(i);
				}
		}

		template<typename T, typename Func>
		inline typename std::enable_if<std::is_scalar<T>::value,T>::type
		omp_coords_range<2>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
#pragma omp parallel for schedule(static) collapse(2) reduction(+: sum)
			for(coord j = r.lower[0]; j < r.upper[0]; j++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
				for(coord k = r.lower[1]; k < r.upper[1]; k++) {
					coords<2> i;
					i[0] = j;
					i[1] = k;
					f(sum, i);
				}
			return sum;
		}

		template<typename T, typename Func>
		inline typename std::enable_if<!std::is_scalar<T>::value,T>::type
		omp_coords_range<2>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
#pragma omp parallel shared(sum)
			{
				T local_sum(zero);
#pragma omp for schedule(static) collapse(2)
				for(coord j = r.lower[0]; j < r.upper[0]; j++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
					for(coord k = r.lower[1]; k < r.upper[1]; k++) {
						coords<2> i;
						i[0] = j;
						i[1] = k;
						f(local_sum, i);
					}
#pragma omp critical
				sum += local_sum;
			}
			return sum;
		}

		template<>
		struct omp_coords_range<3> {
			coords_range<3> r;

			template<typename Func>
			INLINE void for_each(const Func& f) const;

			template<typename T, typename Func>
			INLINE typename std::enable_if<std::is_scalar<T>::value,T>::type
			internal_sum(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			INLINE typename std::enable_if<!std::is_scalar<T>::value,T>::type
			internal_sum(const T& zero, const Func& f) const;
		};

		template<typename Func>
		inline void omp_coords_range<3>::for_each(const Func& f) const {
#pragma omp parallel for schedule(static) collapse(3)
			for(coord j = r.lower[0]; j < r.upper[0]; j++)
				for(coord k = r.lower[1]; k < r.upper[1]; k++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
					for(coord l = r.lower[2]; l < r.upper[2]; l++) {
						coords<3> i;
						i[0] = j;
						i[1] = k;
						i[2] = l;
						f(i);
					}

		}

		template<typename T, typename Func>
		inline typename std::enable_if<std::is_scalar<T>::value,T>::type
		omp_coords_range<3>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
#pragma omp parallel for schedule(static) collapse(3) reduction(+: sum)
			for(coord j = r.lower[0]; j < r.upper[0]; j++)
				for(coord k = r.lower[1]; k < r.upper[1]; k++) 
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
					for(coord l = r.lower[2]; l < r.upper[2]; l++) {
						coords<3> i;
						i[0] = j;
						i[1] = k;
						i[2] = l;
						f(sum, i);
					}
			return sum;
		}

		template<typename T, typename Func>
		inline typename std::enable_if<!std::is_scalar<T>::value,T>::type
		omp_coords_range<3>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
#pragma omp parallel shared(sum)
			{
				T local_sum(zero);
#pragma omp for schedule(static) collapse(3)
				for(coord j = r.lower[0]; j < r.upper[0]; j++)
					for(coord k = r.lower[1]; k < r.upper[1]; k++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
						for(coord l = r.lower[2]; l < r.upper[2]; l++) {
							coords<3> i;
							i[0] = j;
							i[1] = k;
							i[2] = l;
							f(local_sum, i);
						}
#pragma omp critical
				sum += local_sum;
			}
			return sum;
		}
#endif

	}

}

#endif
