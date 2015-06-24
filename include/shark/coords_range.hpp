/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
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

#ifdef SHARK_OMP_SCHED
#ifdef SHARK_RANGE

#define SHARK_RANGE_STEP 256

#endif
#endif

#ifdef __INTEL_COMPILER
#define SHARK_VECTOR_PRAGMA _Pragma(ivdep)
#else 
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
			INLINE typename std::enable_if<d < ndim-1>::type for_eachd_outer_block(const Func& f, coords<ndim>& i) const;
			template<int d, typename Func>
			INLINE typename std::enable_if<d == ndim-1>::type for_eachd_outer_block(const Func& f, coords<ndim>& i) const;
			template<int d, typename Func>
			INLINE typename std::enable_if<d < ndim-1>::type for_eachd_inner_block(const Func& f, coords<ndim>& i, const coords<ndim>& b) const;
			template<int d, typename Func>
			INLINE typename std::enable_if<d == ndim-1>::type for_eachd_inner_block(const Func& f, coords<ndim>& i, const coords<ndim>& b) const;

#ifdef SHARK_RANGE
			template<int d, typename Func>
			INLINE typename std::enable_if<d < ndim-1>::type for_eachr(const Func& f, coords<ndim>& i) const;
			template<int d, typename Func>
			INLINE typename std::enable_if<d == ndim-1>::type for_eachr(const Func& f, coords<ndim>& i) const;
#endif

		public:
			template<typename Func>
			INLINE void for_each(const Func& f) const;
			template<typename Func>
			INLINE void for_each_blocked(const Func& f) const;

#ifdef SHARK_RANGE
			template<typename Func>
			INLINE void for_each_range(const Func& f) const;
#endif

			coords_range<ndim> overlap(coords_range<ndim> other) const;
			bool contains(coords<ndim> i) const;
			bool contains(coords_range<ndim> other) const;
			INLINE coords<ndim> counts() const;
		        INLINE coords<ndim> next(coords<ndim> ii) const;
			INLINE coord count() const;
			INLINE coords<ndim> periodic_equiv(coords<ndim> i) const;
			INLINE coords<ndim> adj(int d, int shift) const;
			coords<ndim+1> stride(coords<ndim> bw = coords<ndim>()) const;
		};

		template<int ndim>
		std::ostream& operator<<(std::ostream& out, const coords_range<ndim>& r);

		// Inline function implementations

#ifdef SHARK_SIMPLE_FOREACH
		template<int ndim> template<typename Func>
		inline void coords_range<ndim>::for_each(const Func& f) const {
                        SHARK_VECTOR_PRAGMA
                        for( auto i = lower; i < upper; i = next(i)) f(i);
		}

#else

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

                        SHARK_VECTOR_PRAGMA
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
#endif

                static const int bs = 64;

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d < ndim-1>::type coords_range<ndim>::for_eachd_outer_block(const Func& f, coords<ndim>& i) const {
			coord ld = lower[d], ud = upper[d];
			for(coord id = ld; id < ud; id+=bs) {
				i[d] = id;
				for_eachd_outer_block<d+1>(f, i);
			}
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d == ndim-1>::type coords_range<ndim>::for_eachd_outer_block(const Func& f, coords<ndim>& i) const {
			coord ld = lower[d], ud = upper[d];
			for(coord id = ld; id < ud; id+=bs) {
				i[d] = id;
                                coords<ndim> b = i;
                                coords<ndim> j = i;
                                for_eachd_inner_block<0>(f, j, b);
			}
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d < ndim-1>::type coords_range<ndim>::for_eachd_inner_block(const Func& f, coords<ndim>& j, const coords<ndim>& b) const {
			coord ld = b[d], ud = std::min(ld+bs,upper[d]);
			for(coord id = ld; id < ud; id++) {
				j[d] = id;
				for_eachd_inner_block<d+1>(f, j, b);
			}
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d == ndim-1>::type coords_range<ndim>::for_eachd_inner_block(const Func& f, coords<ndim>& j, const coords<ndim>& b) const {
			coord ld = b[d], ud = std::min(ld+bs,upper[d]);
                        SHARK_VECTOR_PRAGMA
			for(coord id = ld; id < ud; id++) {
				j[d] = id;
                                f(j);
			}
		}

		template<int ndim> template<typename Func>
		inline void coords_range<ndim>::for_each_blocked(const Func& f) const {
			coords<ndim> i;
			for_eachd_outer_block<0>(f, i);
		}

#ifdef SHARK_RANGE
		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d < ndim-1>::type coords_range<ndim>::for_eachr(const Func& f, coords<ndim>& i) const {
			for(i[d] = lower[d]; i[d] < upper[d]; i[d]++)
				for_eachr<d+1>(f, i);
		}

		template<int ndim> template<int d, typename Func>
		inline typename std::enable_if<d == ndim-1>::type coords_range<ndim>::for_eachr(const Func& f, coords<ndim>& i) const {
			i[d] = lower[d];
			f(i, upper[d] - lower[d]);
		}

		template<int ndim> template<typename Func>
		inline void coords_range<ndim>::for_each_range(const Func& f) const {
			coords<ndim> i;
			for_eachr<0>(f, i);
		}
#endif

		template<int ndim>
		inline coords<ndim> coords_range<ndim>::counts() const {
			return upper - lower;
		}

		template<int ndim>
		inline coords<ndim> coords_range<ndim>::next(coords<ndim> ii) const {
                        ii[ndim-1]++;
		        seq<0,ndim>::for_each([this,&ii](int j) {
                                auto k = ndim - j - 1;
                                if (ii[k] >= upper[k] && k > 0) {
                                    ii[k] = lower[k];
                                    ii[k-1]++;
                                }
                        });
                        return ii;
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

#ifdef SHARK_RANGE
			template<typename Func>
			INLINE void for_each_range(const Func& f) const;

			template<typename T, typename Func>
			INLINE typename std::enable_if<std::is_scalar<T>::value,T>::type
			internal_sum_range(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			INLINE typename std::enable_if<!std::is_scalar<T>::value,T>::type
			internal_sum_range(const T& zero, const Func& f) const;
#endif
		};

		template<typename Func>
		inline void omp_coords_range<1>::for_each(const Func& f) const {
			coord lj = r.lower[0], uj = r.upper[0];
#pragma omp parallel for schedule(static)
                        SHARK_VECTOR_PRAGMA
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
                        SHARK_VECTOR_PRAGMA
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
                                SHARK_VECTOR_PRAGMA
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

#ifdef SHARK_RANGE
	  template<typename Func>
	  inline void omp_coords_range<1>::for_each_range(const Func& f) const {

	    auto start = r.lower[0];
	    auto end = r.upper[0];

#pragma omp parallel for schedule(static)
            SHARK_VECTOR_PRAGMA
	    for(coord j = start; j < end; j+=SHARK_RANGE_STEP) {
	      coords<1> i;
	      i[0] = j;
	      f(i, std::min(j+SHARK_RANGE_STEP,end)-j);
	    }
	  }
	    
	  template<typename T, typename Func>
	  inline typename std::enable_if<std::is_scalar<T>::value,T>::type
	  omp_coords_range<1>::internal_sum_range(const T& zero, const Func& f) const {
	    T sum(zero);
	    
	    auto start = r.lower[0];
	    auto end = r.upper[0];

#pragma omp parallel for schedule(static) reduction(+: sum)
            SHARK_VECTOR_PRAGMA
	    for(coord j = start; j < end; j+=SHARK_RANGE_STEP) {
	      coords<1> i;
	      i[0] = j;
	      f(sum, i, std::min(j+SHARK_RANGE_STEP,end)-j);
	    }
	    return sum;
	  }
	  
	  template<typename T, typename Func>
	  inline typename std::enable_if<!std::is_scalar<T>::value,T>::type
	  omp_coords_range<1>::internal_sum_range(const T& zero, const Func& f) const {
	    T sum(zero);
	    
	    auto start = r.lower[0];
	    auto end = r.upper[0];
	    
#pragma omp parallel shared(sum)
	    {
	      T local_sum(zero);
	      
#pragma omp for schedule(static)
              SHARK_VECTOR_PRAGMA
	      for(coord j = start; j < end; j+=SHARK_RANGE_STEP) {
		coords<1> i;
		i[0] = j;
		f(local_sum, i, std::min(j+SHARK_RANGE_STEP,end)-j);
	      }
#pragma omp critical
	      sum += local_sum;
	    }
	    return sum;
	  }
#endif
	  
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

#ifdef SHARK_RANGE
			template<typename Func>
			INLINE void for_each_range(const Func& f) const;

			template<typename T, typename Func>
			INLINE typename std::enable_if<std::is_scalar<T>::value,T>::type
			internal_sum_range(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			INLINE typename std::enable_if<!std::is_scalar<T>::value,T>::type
			internal_sum_range(const T& zero, const Func& f) const;
#endif
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
                                SHARK_VECTOR_PRAGMA
				for (coord k = lk; k < uk; k++) {
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
                                SHARK_VECTOR_PRAGMA
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
                                        SHARK_VECTOR_PRAGMA
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
		
#ifdef SHARK_RANGE
	  template<typename Func>
	  inline void omp_coords_range<2>::for_each_range(const Func& f) const {
	    
	    auto start = r.lower[1];
	    auto end = r.upper[1];
	    
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(2)
#else
#pragma omp parallel for schedule(static)
#endif
	    for(coord j = r.lower[0]; j < r.upper[0]; j++)
              SHARK_VECTOR_PRAGMA
	      for(coord k = start; k < end; k+=SHARK_RANGE_STEP) {
		coords<2> i;
		i[0] = j;
		i[1] = k;
		f(i, std::min(k+SHARK_RANGE_STEP,end)-k);
	      }
	  }

	  template<typename T, typename Func>
	  inline typename std::enable_if<std::is_scalar<T>::value,T>::type
	  omp_coords_range<2>::internal_sum_range(const T& zero, const Func& f) const {
	    T sum(zero);
	    
	    auto start = r.lower[1];
	    auto end = r.upper[1];
	    
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(2) reduction(+: sum)
#else
#pragma omp parallel for schedule(static) reduction(+: sum)
#endif
	    for(coord j = r.lower[0]; j < r.upper[0]; j++)
              SHARK_VECTOR_PRAGMA
	      for(coord k = start; k < end; k+=SHARK_RANGE_STEP) {
		coords<2> i;
		      i[0] = j;
		      i[1] = k;
		      f(sum, i, std::min(k+SHARK_RANGE_STEP,end)-k);
	      }
	    return sum;
	  }
	  
	  template<typename T, typename Func>
	  inline typename std::enable_if<!std::is_scalar<T>::value,T>::type
	  omp_coords_range<2>::internal_sum_range(const T& zero, const Func& f) const {
	    T sum(zero);
	    
	    auto start = r.lower[1];
	    auto end = r.upper[1];
	    
#pragma omp parallel shared(sum)
	    {
	      T local_sum(zero);
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp for schedule(static) collapse(2)
#else
#pragma omp for schedule(static)
#endif
	      for(coord j = r.lower[0]; j < r.upper[0]; j++)
                SHARK_VECTOR_PRAGMA
		for(coord k = start; k < end; k+=SHARK_RANGE_STEP) {
		  coords<2> i;
		  i[0] = j;
		  i[1] = k;
		  f(local_sum, i, std::min(k+SHARK_RANGE_STEP,end)-k);
		}
#pragma omp critical
	      sum += local_sum;
	    }
	    return sum;
	  }
#endif
	  
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

#ifdef SHARK_RANGE
			template<typename Func>
			INLINE void for_each_range(const Func& f) const;

			template<typename T, typename Func>
			INLINE typename std::enable_if<std::is_scalar<T>::value,T>::type
			internal_sum_range(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			INLINE typename std::enable_if<!std::is_scalar<T>::value,T>::type
			internal_sum_range(const T& zero, const Func& f) const;
#endif
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
                                        SHARK_VECTOR_PRAGMA
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
                                        SHARK_VECTOR_PRAGMA
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
                                                SHARK_VECTOR_PRAGMA
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
		
#ifdef SHARK_RANGE
	  template<typename Func>
	  inline void omp_coords_range<3>::for_each_range(const Func& f) const {
	    
	    auto start = r.lower[2];
	    auto end = r.upper[2];
	    
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(3)
#else
#pragma omp parallel for schedule(static)
#endif
	    for(coord j = r.lower[0]; j < r.upper[0]; j++)
	      for(coord k = r.lower[1]; k < r.upper[1]; k++)
                SHARK_VECTOR_PRAGMA
		for(coord l = start; l < end; l+=SHARK_RANGE_STEP) {
		  coords<3> i;
		  i[0] = j;
		  i[1] = k;
		  i[2] = l;
		  f(i, std::min(l+SHARK_RANGE_STEP,end)-l);
		}
	    
	  }
	  
	  template<typename T, typename Func>
	  inline typename std::enable_if<std::is_scalar<T>::value,T>::type
	  omp_coords_range<3>::internal_sum_range(const T& zero, const Func& f) const {
	    T sum(zero);
	    
	    auto start = r.lower[2];
	    auto end = r.upper[2];
	    
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp parallel for schedule(static) collapse(3) reduction(+: sum)
#else
#pragma omp parallel for schedule(static) reduction(+: sum)
#endif
	    for(coord j = r.lower[0]; j < r.upper[0]; j++)
	      for(coord k = r.lower[1]; k < r.upper[1]; k++) 
                SHARK_VECTOR_PRAGMA
		for(coord l = start; l < end; l+=SHARK_RANGE_STEP) {
		  coords<3> i;
		  i[0] = j;
		  i[1] = k;
		  i[2] = l;
		  f(sum, i, std::min(l+SHARK_RANGE_STEP,end)-l);
		}
	    return sum;
	  }
	  
	  template<typename T, typename Func>
	  inline typename std::enable_if<!std::is_scalar<T>::value,T>::type
	  omp_coords_range<3>::internal_sum_range(const T& zero, const Func& f) const {
	    T sum(zero);
	    
	    auto start = r.lower[2];
	    auto end = r.upper[2];
	    
#pragma omp parallel shared(sum)
	    {
	      T local_sum(zero);
#ifdef SHARK_THREAD_BLOCK_DIST
#pragma omp for schedule(static) collapse(3)
#else
#pragma omp for schedule(static) 
#endif
	      for(coord j = r.lower[0]; j < r.upper[0]; j++)
		for(coord k = r.lower[1]; k < r.upper[1]; k++)
                  SHARK_VECTOR_PRAGMA
		  for(coord l = start; l < end; l+=SHARK_RANGE_STEP) {
		    coords<3> i;
		    i[0] = j;
		    i[1] = k;
		    i[2] = l;
		    f(local_sum, i, std::min(l+SHARK_RANGE_STEP,end)-l);
		  }
#pragma omp critical
	      sum += local_sum;
	    }
	    return sum;
	  }
#endif

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
