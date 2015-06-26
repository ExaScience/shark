/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_DOMAIN_HPP
#define __SHARK_DOMAIN_HPP

#include <array>                       // std::array
#include <bitset>                      // std::bitset
#include <memory>                      // std::unique_ptr
#include <vector>                      // std::vector
#include <ostream>                     // std::ostream
#if defined(SHARK_TBB_SCHED)
#include <tbb/tbb_stddef.h>            // tbb::split
#include <tbb/partitioner.h>           // tbb::affinity_partitioner
#include <tbb/parallel_for.h>          // tbb::parallel_for
#include <tbb/parallel_reduce.h>       // tbb::parallel_reduce
#if !defined(SHARK_THREAD_BLOCK_DIST)
#include <tbb/blocked_range.h>         // tbb::blocked_range
#endif
#endif

#include "common.hpp"
#include "coords.hpp"
#include "coords_range.hpp"
#include "group.hpp"
#include "future.hpp"

namespace shark {

	namespace ndim {

		template<int ndim>
		class Domain {
		public:
			typedef std::array<int,ndim> pcoords;
			typedef std::array<std::vector<coord>,ndim> dists;

			/**
			 * The group of processes that are hosting this domain.
			 */
			const Group& group;

			/**
			 * The total number of elements for every dimension.
			 */
			const coords<ndim> n;

			/**
			 * The number of processes for every dimension.
			 */
			const pcoords np;

			/**
			 * The distribution of the elements for every dimension.
			 */
			const dists nd;

			/**
			 * Construct a new domain (local).
			 * The data will be evenly distributed across all processes.
			 * @param n the number of elements for every dimension
			 * @param np the number of processes for every dimension (0 = auto)
			 * \verbatim nprocs = prod{d=0..ndim-1} np[d] \endverbatim
			 */
			Domain(const Group& group, coords<ndim> n, pcoords np = pcoords());

			/**
			 * Construct a new domain (local).
			 * The given distribution will be used
			 * @param n the number of elements for every dimension
			 * @param np the number of processes for every dimension.
			 * \verbatim nprocs = prod{d=0..ndim-1} np[d] \endverbatim
			 * @param dist specifyies the distribution for each dimension
			 * \verbatim for each d=0..ndim-1, dist[d][0] = 0, dist[d][np[d]] = n[d], dist[d][i] <= dist[d][j] when i < j \endverbatim
			 */
			Domain(const Group& group, coords<ndim> n, pcoords np, dists nd);

			/**
			 * Destruct a domain (local).
			 */
			~Domain();

			Domain(const Domain&) = delete;
			Domain& operator=(const Domain&) = delete;

		private:
			const std::array<int,ndim+1> b;

#if defined(SHARK_PTHREAD_SCHED) || defined(SHARK_OMP_SCHED) && defined(SHARK_OMP_TDIST)
			const std::vector<coords_range<ndim>> tdist;
			std::vector<coords_range<ndim>> tdistribution() const;
#elif defined(SHARK_TBB_SCHED)
			const std::unique_ptr<tbb::affinity_partitioner> ap;
			const int nwork;
#endif
			static void adjustProcs(int nprocs, pcoords& np);
			static std::array<int,ndim+1> base(pcoords np);
			static dists distribution(coords<ndim> n, pcoords np);
			bool consistentDistribution() const;

			static void findi(const std::vector<coord>& dist, coord i, int& id, coord& off);

		public:
			/**
			 * Get process coordinates for a process (local).
			 * @param id the process identifier (default: group.procid)
			 */
			INLINE pcoords indexp(int id) const;
			INLINE pcoords indexp() const;

			/**
			 * Get process identifier for process coordinates (local).
			 * @param pc the process coordinates
			 */
			INLINE int pindex(pcoords ip) const; 

			/**
			 * Get process identifier after a process shift (local).
			 * @param d the dimension to shift
			 * @param disp the amount to shift with
			 * @param pd wrap around periodically (else returns MPI_PROC_NULL or -1)
			 * @param id the process identifier to start from (default: group.procid)
			 */
			int shiftd(int d, int disp, bool pd, int id) const;
			INLINE int shiftd(int d, int disp, bool pd) const;

			/**
			 * Get coordinates of the array region owned by a process (local).
			 * @param id the process identifier (default: group.procid)
			 */
			coords_range<ndim> local(int id) const;
			INLINE coords_range<ndim> local() const;

			/**
			 * Get coordinates of the total domain region (local).
			 */
			INLINE coords_range<ndim> total() const;

			/**
			 * Check whether a process owns any data (local).
			 * @param id the process identifier (default: group.procid)
			 */
			INLINE bool hasData(int id) const;
			INLINE bool hasData() const;

			/**
			 * Return the number of cells per dimension owned by a process (local).
			 * @param id the process identifier (default: group.procid)
			 */
			INLINE coords<ndim> count(int id) const;
			INLINE coords<ndim> count() const;

			/**
			 * Output the distribution of this global array (local).
			 * @param out the output stream to write to
			 */
			void outputDistribution(std::ostream& out) const;

			/**
			 * Efficiently find which process owns a given point (local).
			 * @param i the coordinates of the point to find
			 * @param ip the process coordinates (output)
			 * @param off the offset of the point at the process (output)
			 */
			void find(coords<ndim> i, pcoords& ip, coords<ndim>& off) const;

			/**
			 * Check whether this domain is equal to another one (local).
			 * @param other the domain to compare to
			 */
			bool operator==(const Domain<ndim>& other) const;

			/**
			 * Check whether this domain is different from another one (local).
			 * @param other the domain to compare to
			 */
			INLINE bool operator!=(const Domain<ndim>& other) const;

			/**
			 * Check whether this domain is equivalent to another one (local).
			 * Equivalent means that it distributes the same number of elements
			 * over the same group of processes.
			 * @param other the domain to compare to
			 */
			bool equiv(const Domain<ndim>& other) const;

			/**
			 * Synchronize the processes of this domain (collective).
			 */
			void sync() const;

			/**
			 * Apply elemental function onto the elements of range
			 * @param r the range where function is applied (default: total())
			 * @param f the function to apply
			 */
			template<typename Func>
			INLINE void for_each(const Func& f) const;
			template<typename Func>
			void for_each(coords_range<ndim> r, const Func& f) const;
#ifdef SHARK_RANGE
                        template<typename Func>
                        INLINE void for_each_range(const Func& f) const;
                        template<typename Func>
                        void for_each_range(coords_range<ndim> r, const Func& f) const;
#endif

			/**
			 * Sum elemental function over the elements of range
			 */
			template<typename T, typename Func>
			INLINE T sum(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			T sum(coords_range<ndim> r, const T& zero, const Func& f) const;

#ifdef SHARK_RANGE
                        template<typename T, typename Func>
                        INLINE T sum_range(const T& zero, const Func& f) const;
                        template<typename T, typename Func>
                        T sum_range(coords_range<ndim> r, const T& zero, const Func& f) const;
#endif

			template<typename T, typename Func>
			INLINE Future<T> isum(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			Future<T> isum(coords_range<ndim> r, const T& zero, const Func& f) const;

			template<typename T, typename Func>
			INLINE T internal_sum(const T& zero, const Func& f) const;
			template<typename T, typename Func>
			T internal_sum(coords_range<ndim> r, const T& zero, const Func& f) const;

#ifdef SHARK_RANGE
                        template<typename T, typename Func>
                        INLINE T internal_sum_range(const T& zero, const Func& f) const;
                        template<typename T, typename Func>
                        T internal_sum_range(coords_range<ndim> r, const T& zero, const Func& f) const;
#endif

			class ProcessOverlap {
				const Domain<ndim>& dom;
				const coords_range<ndim> range;

				pcoords lip, uip;
				coords<ndim> loff, uoff;

				template<int d,typename Op>
				INLINE typename std::enable_if<d < ndim>::type visitd(const Op& op, pcoords& ip, coords_range<ndim>& i);

				template<int d,typename Op>
				INLINE typename std::enable_if<d == ndim>::type visitd(const Op& op, pcoords& ip, coords_range<ndim>& i);

			public:
				ProcessOverlap(const Domain<ndim>& dom, coords_range<ndim> range);

				template<typename Op>
				void visit(const Op& op);
			};

		};

		// Inline function implementations
		
		template<int ndim>
		inline typename Domain<ndim>::pcoords Domain<ndim>::indexp(int id) const
		{
			pcoords ip;

			seq<0,ndim>::for_each([&id,&ip,this](int d)
			{
				ip[d] = id / b[d+1];
				id = id % b[d+1];
			}
			);

			return ip;
		}

		template<int ndim>
		inline typename Domain<ndim>::pcoords Domain<ndim>::indexp() const {
			return indexp(group.procid);
		}

		template<int ndim>
		inline int Domain<ndim>::pindex(pcoords ip) const {
			int id = 0;

			seq<0,ndim>::for_each([&id,&ip,this](int d)
			{
				id += ip[d] * b[d+1];
			});

			return id;
		}
		
		template<int ndim>
		inline int Domain<ndim>::shiftd(int d, int disp, bool pd) const {
			return shiftd(d, disp, pd, group.procid);
		}

		template<int ndim>
		inline coords_range<ndim> Domain<ndim>::local() const {
			return local(group.procid);
		}


		template<int ndim>
		inline coords_range<ndim> Domain<ndim>::total() const {
			const coords_range<ndim> r = {{{}}, n};
			return r;
		}
		
		template<int ndim>
		inline bool Domain<ndim>::hasData(int id) const {
			const coords_range<ndim> r = local(id);
			return seq<0,ndim>::all_of([&r](int d) {
				return r.lower[d] < r.upper[d];
			});
		}

		template<int ndim>
		inline bool Domain<ndim>::hasData() const {
			return hasData(group.procid);
		}

		template<int ndim>
		inline coords<ndim> Domain<ndim>::count(int id) const {
			const coords_range<ndim> r = local(id);
			return r.upper - r.lower;
		}

		template<int ndim>
		inline coords<ndim> Domain<ndim>::count() const {
			return count(group.procid);
		}

		template<int ndim>
		inline bool Domain<ndim>::operator!=(const Domain<ndim>& other) const {
			return !(*this == other);
		}

		// Generic Domain members
		
		template<int ndim> template<typename Func>
		inline void Domain<ndim>::for_each(const Func& f) const {
			for_each(total(), f);
		}
		
		template<int ndim> template<typename Func>
		void Domain<ndim>::for_each(coords_range<ndim> r, const Func& f) const {
#if defined(SHARK_SER_SCHED)
			local().overlap(r).for_each(f);
#elif defined(SHARK_PTHREAD_SCHED)
			ThreadWork([this,&f,r](int k) {
				tdist[k].overlap(r).for_each(f);
			});
#elif defined(SHARK_OMP_SCHED)
#if defined(SHARK_OMP_TDIST)
#pragma omp parallel for schedule(static)
			for(std::size_t k = 0; k < tdist.size(); k++)
				tdist[k].overlap(r).for_each(f);
#else
			omp_coords_range<ndim> omp;
			omp.r = local().overlap(r);
			omp.for_each(f);
#endif
#elif defined(SHARK_TBB_SCHED)
			r = local().overlap(r);
#if defined(SHARK_THREAD_BLOCK_DIST)
			split_range<ndim> sr(local(), r.count()/nwork);
			tbb::parallel_for(sr, [&r,&f](const split_range<ndim>& sr) {
				sr.range().overlap(r).for_each(f);
			}, *ap);
#else
			tbb::blocked_range<coord> br(local().lower[0], local().upper[0]);
			tbb::parallel_for(br, [&r,&f](const tbb::blocked_range<coord>& br) {
				coords_range<ndim> lr = r;
				if(lr.lower[0] < br.begin())
					lr.lower[0] = br.begin();
				if(lr.upper[0] > br.end())
					lr.upper[0] = br.end();
				lr.for_each_blocked(f);
				//lr.for_each(f);
			}, *ap);
#endif
#else
#error "No scheduler for_each"
#endif
		}

#ifdef SHARK_RANGE
		template<int ndim> template<typename Func>
		inline void Domain<ndim>::for_each_range(const Func& f) const {
			for_each_range(total(), f);
		}
		
		template<int ndim> template<typename Func>
		void Domain<ndim>::for_each_range(coords_range<ndim> r, const Func& f) const {
#if defined(SHARK_SER_SCHED)
			local().overlap(r).for_each_range(f);
#elif defined(SHARK_PTHREAD_SCHED)
			ThreadWork([this,&f,r](int k) {
				tdist[k].overlap(r).for_each_range(f);
			});
#elif defined(SHARK_OMP_SCHED)
			omp_coords_range<ndim> omp;
			omp.r = local().overlap(r);
			omp.for_each_range(f);
#else
#error "No scheduler for_each_range"
#endif
		}
#endif

		template<int ndim> template<typename T, typename Func>
		inline T Domain<ndim>::internal_sum(const T& zero, const Func& f) const {
			return internal_sum(total(), zero, f);
		}
		
		template<int ndim> template<typename T, typename Func>
		T Domain<ndim>::internal_sum(coords_range<ndim> r, const T& zero, const Func& f) const
		{
#if defined(SHARK_SER_SCHED)
			T sum(zero);
			local().overlap(r).for_each([&f,&sum](coords<ndim> i) {
				f(sum, i);
			});
			return sum;
#elif defined(SHARK_PTHREAD_SCHED)
			std::unique_ptr<T>* tsum = new std::unique_ptr<T>[tdist.size()];
			ThreadWork([this,&f,&tsum,&zero,r](int k) {
				tsum[k].reset(new T(zero));
				tdist[k].overlap(r).for_each([&f,&tsum,k](coords<ndim> i) {
					f(*tsum[k], i);
				});
			});
			std::unique_ptr<T> sum(std::move(tsum[0]));
			for(std::size_t k=1; k < tdist.size(); k++)
				*sum += *tsum[k];
			delete[] tsum;
			return *sum;
#elif defined(SHARK_OMP_SCHED)
#if defined(SHARK_OMP_TDIST)
			std::unique_ptr<T>* tsum = new std::unique_ptr<T>[tdist.size()];
#pragma omp parallel for schedule(static)
			for(std::size_t k = 0; k < tdist.size(); k++) {
				tsum[k].reset(new T(zero));
				tdist[k].overlap(r).for_each([&f,&tsum,k](coords<ndim> i) {
					f(*tsum[k], i);
				});
			}
			std::unique_ptr<T> sum(std::move(tsum[0]));
			for(std::size_t k=1; k < tdist.size(); k++)
				*sum += *tsum[k];
			delete[] tsum;
			return *sum;
#else
			omp_coords_range<ndim> omp;
			omp.r = local().overlap(r);
			return omp.internal_sum(zero, f);
#endif
#elif defined(SHARK_TBB_SCHED)
			r = local().overlap(r);
#if defined(SHARK_THREAD_BLOCK_DIST)
			struct Adder {
				const Func& f;
				const T& zero;
				const coords_range<ndim>& r;
				T sum;
				Adder(const Func& f, const T& zero, const coords_range<ndim>& r) : f(f), zero(zero), r(r), sum(zero) {}
				Adder(Adder& adder, tbb::split) : f(adder.f), zero(adder.zero), r(adder.r), sum(zero) {}
				void operator()(const split_range<ndim>& sr) {
					sr.range().overlap(r).for_each([this](coords<ndim> i) {
						f(sum, i);
					});
				}
				void join(Adder& rhs) { sum += rhs.sum; }
			} adder(f, zero, r);
			split_range<ndim> sr(local(), r.count()/nwork);
			tbb::parallel_reduce(sr, adder, *ap);
			return adder.sum;
#else
			struct Adder {
				const Func& f;
				const T& zero;
				const coords_range<ndim>& r;
				T sum;
				Adder(const Func& f, const T& zero, const coords_range<ndim>& r) 
					: f(f), zero(zero), r(r), sum(zero) {}
				Adder(Adder& adder, tbb::split)
					: f(adder.f), zero(adder.zero), r(adder.r), sum(zero) {}

				void operator()(const tbb::blocked_range<coord>& br) {
					coords_range<ndim> lr = r;
					if(lr.lower[0] < br.begin())
						lr.lower[0] = br.begin();
					if(lr.upper[0] > br.end())
						lr.upper[0] = br.end();
					lr.for_each([this](coords<ndim> i) {
						f(sum, i);
					});
				}
				void join(Adder& rhs) { sum += rhs.sum; }
			} adder(f, zero, r);
			tbb::blocked_range<coord> br(local().lower[0], local().upper[0]);
			tbb::parallel_reduce(br, adder, *ap);
			return adder.sum;
#endif
#else
#error "No scheduler internal_sum"
#endif
		}

#ifdef SHARK_RANGE
		template<int ndim> template<typename T, typename Func>
		inline T Domain<ndim>::internal_sum_range(const T& zero, const Func& f) const {
			return internal_sum_range(total(), zero, f);
		}
		
		template<int ndim> template<typename T, typename Func>
		T Domain<ndim>::internal_sum_range(coords_range<ndim> r, const T& zero, const Func& f) const {
#if defined(SHARK_SER_SCHED)
			T sum(zero);
			local().overlap(r).for_each_range([&f,&sum](coords<ndim> i, coord len) {
				f(sum, i, len);
			});
			return sum;
#elif defined(SHARK_PTHREAD_SCHED)
			std::vector<T> tsum(nthrds, zero);
			ThreadWork([this,&f,&tsum,r](int k) {
				tdist[k].overlap(r).for_each_range([&f,&tsum,k](coords<ndim> i, coord len) {
					f(tsum[k], i, len);
				});
			});
			for(int k=1; k < nthrds; k++)
				tsum[0] += tsum[k];
			return tsum[0];
#elif defined(SHARK_OMP_SCHED)
			omp_coords_range<ndim> omp;
			omp.r = local().overlap(r);
			return omp.internal_sum_range(zero, f);
#else
#error "No scheduler internal_sum_range"
#endif
		}
#endif

		template<int ndim> template<typename T, typename Func>
		inline T Domain<ndim>::sum(const T& zero, const Func& f) const {
			return sum(total(), zero, f);
		}

		template<int ndim> template<typename T, typename Func>
		T Domain<ndim>::sum(coords_range<ndim> r, const T& zero, const Func& f) const
		{
			return group.external_sum(internal_sum(r, zero, f));
		}

		template<int ndim> template<typename T, typename Func>
		inline Future<T> Domain<ndim>::isum(const T& zero, const Func& f) const {
			return isum(total(), zero, f);
		}

		template<int ndim> template<typename T, typename Func>
		Future<T> Domain<ndim>::isum(coords_range<ndim> r, const T& zero, const Func& f) const {
			return group.external_isum(internal_sum(r, zero, f));
		}
  
#ifdef SHARK_RANGE
                template<int ndim> template<typename T, typename Func>
                inline T Domain<ndim>::sum_range(const T& zero, const Func& f) const {
                  return sum_range(total(), zero, f);
                }
                
                template<int ndim> template<typename T, typename Func>
                T Domain<ndim>::sum_range(coords_range<ndim> r, const T& zero, const Func& f) const {
                  return group.external_sum(internal_sum_range(r, zero, f));
                }
#endif

		template<int ndim> template<typename Op>
		void Domain<ndim>::ProcessOverlap::visit(const Op& op) {
			pcoords ip;
			coords_range<ndim> i;
			visitd<0>(op, ip, i);
		}

		template<int ndim> template<int d,typename Op>
		inline typename std::enable_if<d < ndim>::type Domain<ndim>::ProcessOverlap::visitd(const Op& op, pcoords& ip, coords_range<ndim>& i) {
			for(ip[d] = lip[d]; ip[d] <= uip[d]; ip[d]++) {
				i.lower[d] = ip[d] == lip[d] ? range.lower[d] : dom.nd[d][ip[d]];
				i.upper[d] = ip[d] == uip[d] ? range.upper[d] : dom.nd[d][ip[d]+1];
				visitd<d+1>(op, ip, i);
			}
		}

		template<int ndim> template<int d,typename Op>
		inline typename std::enable_if<d == ndim>::type Domain<ndim>::ProcessOverlap::visitd(const Op& op, pcoords& ip, coords_range<ndim>& i) {
			// Target
			int id = dom.pindex(ip);
			op(id, i);
		}
	}
}

#endif
