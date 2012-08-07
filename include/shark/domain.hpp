#ifndef __SHARK_DOMAIN_HPP
#define __SHARK_DOMAIN_HPP

#include <array>                       // std::array
#include <bitset>                      // std::bitset
#include <memory>                      // std::unique_ptr
#include <vector>                      // std::vector
#include <ostream>                     // std::ostream
#include "common.hpp"
#include "coords.hpp"
#include "coords_range.hpp"
#include "group.hpp"

namespace shark {

	namespace ndim {

		template<int ndim>
		class Domain {
			friend class Region<ndim>;

		public:
			typedef std::array<int,ndim> pcoords;
			typedef std::array<std::vector<coord>,ndim> dists;
			typedef std::bitset<ndim> periods;

			/**
			 * The group of processes that are hosting this domain.
			 */
			const Group& group;

			/**
			 * The total number of elements for every dimension.
			 */
			const coords<ndim> n;

			/**
			 * The periodicity of the domain for every dimension.
			 */
			const periods pd;

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
			 * @param pd the periodicity of the domain for every dimension
			 * @param np the number of processes for every dimension (0 = auto)
			 * \verbatim nprocs = prod{d=0..ndim-1} np[d] \endverbatim
			 */
			Domain(const Group& group, coords<ndim> n, periods pd = periods(), pcoords np = pcoords());

			/**
			 * Construct a new domain (local).
			 * The given distribution will be used
			 * @param n the number of elements for every dimension
			 * @param np the number of processes for every dimension.
			 * \verbatim nprocs = prod{d=0..ndim-1} np[d] \endverbatim
			 * @param dist specifyies the distribution for each dimension
			 * \verbatim for each d=0..ndim-1, dist[d][0] = 0, dist[d][np[d]] = n[d], dist[d][i] <= dist[d][j] when i < j \endverbatim
			 */
			Domain(const Group& group, coords<ndim> n, periods pd, pcoords np, dists nd);

			/**
			 * Destruct a domain (local).
			 */
			~Domain();

			Domain(const Domain&) = delete;
			Domain& operator=(const Domain&) = delete;

		private:
			const std::array<int,ndim+1> b;
#if defined(SHARK_PTHREAD_SCHED)
			const std::vector<coords_range<ndim>> tdist;
			std::vector<coords_range<ndim>> tdistribution() const;
#elif defined(SHARK_TBB_SCHED)
			const std::unique_ptr<tbb::affinity_partitioner> ap;
#endif
			static void adjustProcs(int nprocs, pcoords& np);
			static std::array<int,ndim+1> base(pcoords np);
			static dists distribution(coords<ndim> n, pcoords np);
			bool consistentDistribution() const;

			static void findi(const std::vector<coord>& dist, coord i, int& id, coord& off);

			template <typename T>
			T external_add_reduce(T&& val) const;

#ifdef SHARK_ASYNC
			template <typename T>
			Future<T> external_iadd_reduce(T&& val) const;
#endif

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
			 * Get process identifier after a process shift (local). This will wrap
			 * around in the period case or return MPI_PROC_NULL in the non-period case.
			 * @param d the dimension to shift
			 * @param disp the amount to shift with
			 * @param id the process identifier to start from (default: group.procid)
			 */
			int shiftd(int d, int disp, int id) const;
			INLINE int shiftd(int d, int disp) const;

			/**
			 * Get coordinates of the array region owned by a process (local).
			 * @param id the process identifier (default: group.procid)
			 */
			INLINE coords_range<ndim> local(int id) const;
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
			 * with the same periodicity over the same group of processes.
			 * @param other the domain to compare to
			 */
			bool equiv(const Domain<ndim>& other) const;

			/**
			 * Synchronize the processes of this domain (collective).
			 */
			void sync() const;

			/**
			 * Apply elemental function onto the elements of the domain
			 */
			template<typename Func>
			void for_each(const Func& f) const;

			/**
			 * Apply elemental function onto the elements of range
			 */
			template<typename Func>
			void for_each(coords_range<ndim> r, const Func& f) const;
		};

		// Inline function implementations
		
		template<int ndim>
		inline typename Domain<ndim>::pcoords Domain<ndim>::indexp(int id) const {
			pcoords ip;
			seq<0,ndim>::for_each([&id,&ip,this](int d){
				ip[d] = id / b[d+1];
				id = id % b[d+1];
			});
			return ip;
		}

		template<int ndim>
		inline typename Domain<ndim>::pcoords Domain<ndim>::indexp() const {
			return indexp(group.procid);
		}

		template<int ndim>
		inline int Domain<ndim>::pindex(pcoords ip) const {
			int id = 0;
			seq<0,ndim>::for_each([&id,&ip,this](int d) {
				id += ip[d] * b[d+1];
			});
			return id;
		}
		
		template<int ndim>
		inline int Domain<ndim>::shiftd(int d, int disp) const {
			return shiftd(d, disp, group.procid);
		}

		template<int ndim>
		inline coords_range<ndim> Domain<ndim>::local(int id) const {
			pcoords ip = indexp(id);
			coords_range<ndim> r;
			seq<0,ndim>::for_each([this,&r,&ip](int d) {
				r.lower[d] = nd[d][ip[d]];
				r.upper[d] = nd[d][ip[d]+1];
			});
			return r;
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
		void Domain<ndim>::for_each(const Func& f) const {
			for_each(total(), f);
		}

		template<int ndim> template<typename Func>
		void Domain<ndim>::for_each(coords_range<ndim> r, const Func& f) const {
			r = local().overlap(r);
#if defined(SHARK_SER_SCHED)
			r.for_each(f);
#else
#error "No scheduler for_each"
#endif
		}
	}
}

#endif
