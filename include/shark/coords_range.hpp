/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_COORDS_RANGE_HPP
#define __SHARK_COORDS_RANGE_HPP

#include <ostream>                     // std::ostream
#if defined(SHARK_TBB_SCHED)
#include <tbb/tbb_stddef.h>            // tbb::split
#endif

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

		  template<int d, typename Func>
		  INLINE typename std::enable_if<d < ndim-1>::type for_eachr(const Func& f, const coords<ndim>& i) const;

		  template<int d, typename Func>
		  INLINE typename std::enable_if<d == ndim-1>::type for_eachr(const Func& f, const coords<ndim>& i) const;

		public:
			template<typename Func>
			INLINE void for_each(const Func& f) const;

		  template<typename Func>
		  INLINE void for_each_range(const Func& f) const;

			coords_range<ndim> overlap(coords_range<ndim> other) const;
			bool contains(coords<ndim> i) const;
			bool contains(coords_range<ndim> other) const;
			INLINE coords<ndim> counts() const;
			INLINE coord count() const;
			INLINE coords<ndim> periodic_equiv(coords<ndim> i) const;
			INLINE coords<ndim> adj(int d, int shift) const;
			coords<ndim+1> stride(coords<ndim> bw = coords<ndim>()) const;
		};

		template<int ndim>
		std::ostream& operator<<(std::ostream& out, const coords_range<ndim>& r);

		// Inline function implementations

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d < ndim-1>::type coords_range<ndim>::for_eachd(const Func& f, coords<ndim>& i) const {
			coord ld = lower[d], ud = upper[d];
			for(coord id = ld; id < ud; id++) {
				i[d] = id;
				for_eachd<d+1>(f, i);
			}
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d == ndim-1>::type coords_range<ndim>::for_eachd(const Func& f, coords<ndim>& i) const {
			coord ld = lower[d], ud = upper[d];
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
			for(coord id = ld; id < ud; id++) {
				i[d] = id;
				f(i);
			}
		}

		template<int ndim> template<typename Func>
		inline void coords_range<ndim>::for_each(const Func& f) const {
			coords<ndim> i;
			for_eachd<0>(f, i);
		}

	  template<int ndim> template<int d, typename Func>
	  inline typename std::enable_if<d < ndim-1>::type coords_range<ndim>::for_eachr(const Func& f, const coords<ndim>& i) const {
	    coords<ndim> ic(i);
	    const auto low = lower[d];
	    const auto up = upper[d];
	    for (auto j = low; j < up; j++) {
	      ic[d] = j;
	      for_eachr<d+1>(f, ic);
	    }
	  }
	    
	    template<int ndim> template<int d, typename Func>
	    inline typename std::enable_if<d == ndim-1>::type coords_range<ndim>::for_eachr(const Func& f, const coords<ndim>& i) const {
	      f(i, upper[d] - lower[d]);
	    }
	  
	  template<int ndim> template<typename Func>
	  inline void coords_range<ndim>::for_each_range(const Func& f) const {
	    coords<ndim> i(lower);
	    for_eachr<0>(f, i);
	  }

		template<int ndim>
		inline coords<ndim> coords_range<ndim>::counts() const {
			return upper - lower;
		}

		template<int ndim>
		inline coord coords_range<ndim>::count() const {
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

#if defined(SHARK_OMP_SCHED) && !defined(SHARK_OMP_TDIST)

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
			coord lj = r.lower[0], uj = r.upper[0];
#pragma omp parallel for schedule(static)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
			for(coord j = lj; j < uj; j++) {
				coords<1> i;
				i[0] = j;
				f(i);
			}
		}

		template<typename T, typename Func>
		inline typename std::enable_if<std::is_scalar<T>::value,T>::type
		omp_coords_range<1>::internal_sum(const T& zero, const Func& f) const {
			T sum(zero);
			coord lj = r.lower[0], uj = r.upper[0];
#pragma omp parallel for schedule(static) reduction(+: sum)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
			for(coord j = lj; j < uj; j++) {
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
				coord lj = r.lower[0], uj = r.upper[0];
#pragma omp for schedule(static)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
				for(coord j = lj; j < uj; j++) {
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
			coord lj = r.lower[0], uj = r.upper[0], lk = r.lower[1], uk = r.upper[1];
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(2)
#else
#pragma omp parallel for schedule(static)
#endif
			for(coord j = lj; j < uj; j++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
				for(coord k = lk; k < uk; k++) {
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
			coord lj = r.lower[0], uj = r.upper[0], lk = r.lower[1], uk = r.upper[1];
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(2) reduction(+: sum)
#else
#pragma omp parallel for schedule(static) reduction(+: sum)
#endif
			for(coord j = lj; j < uj; j++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
				for(coord k = lk; k < uk; k++) {
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
				coord lj = r.lower[0], uj = r.upper[0], lk = r.lower[1], uk = r.upper[1];
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp for schedule(static) collapse(2)
#else
#pragma omp for schedule(static)
#endif
				for(coord j = lj; j < uj; j++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
					for(coord k = lk; k < uk; k++) {
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
			coord lj = r.lower[0], uj = r.upper[0], lk = r.lower[1], uk = r.upper[1], ll = r.lower[2], ul = r.upper[2];
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(3)
#else
#pragma omp parallel for schedule(static)
#endif
			for(coord j = lj; j < uj; j++)
				for(coord k = lk; k < uk; k++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
					for(coord l = ll; l < ul; l++) {
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
			coord lj = r.lower[0], uj = r.upper[0], lk = r.lower[1], uk = r.upper[1], ll = r.lower[2], ul = r.upper[2];
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(3) reduction(+: sum)
#else
#pragma omp parallel for schedule(static) reduction(+: sum)
#endif
			for(coord j = lj; j < uj; j++)
				for(coord k = lk; k < uk; k++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
					for(coord l = ll; l < ul; l++) {
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
				coord lj = r.lower[0], uj = r.upper[0], lk = r.lower[1], uk = r.upper[1], ll = r.lower[2], ul = r.upper[2];
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp for schedule(static) collapse(3)
#else
#pragma omp for schedule(static) 
#endif
				for(coord j = lj; j < uj; j++)
					for(coord k = lk; k < uk; k++)
#ifdef SHARK_VECTOR_PRAGMA
#pragma ivdep
#endif
						for(coord l = ll; l < ul; l++) {
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
		
#elif defined(SHARK_TBB_SCHED) && defined(SHARK_THREAD_BLOCK_DIST)

		template<int ndim>
		class split_range {
		private:
			coords_range<ndim> r;
			coord grainsize;
		public:
			split_range(coords_range<ndim> r, coord grainsize);
			split_range(const split_range<ndim>& other) = default;
			split_range(split_range<ndim>& other, tbb::split);
			~split_range() = default;
			void split(split_range<ndim>& left, split_range<ndim>& right) const;
			bool is_divisible() const;
			bool empty() const;
			INLINE const coords_range<ndim>& range() const;
		};

		template<int ndim>
		inline const coords_range<ndim>& split_range<ndim>::range() const {
			return r;
		}

#endif

	}

}

#endif
