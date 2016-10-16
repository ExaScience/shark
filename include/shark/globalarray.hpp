/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_GLOBALARRAY_HPP
#define __SHARK_GLOBALARRAY_HPP

#include <array>                       // std::array
#include <cstddef>                     // std::size_t
#include <memory>                      // std::unique_ptr
#include <cassert>                     // assert
#include "common.hpp"
#include "boundary.hpp"
#include "coords.hpp"
#include "coords_range.hpp"
#include "future.hpp"
#include "globals.hpp"
//#include <mpi.h>

namespace shark {

	namespace ndim {

		template<int,typename>
		class GAImpl;

		template<int,typename>
		class GADest;

		template<int,typename>
		class GARef;

		template<int,typename>
		class GABuf;

		template<int ndim,typename T>
                void get(GlobalArray<ndim, T> &tgt, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);
		template<int ndim,typename T>
                void put(GlobalArray<ndim, T> &tgt, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);
		template<int ndim,typename T>
                void accumulate(GlobalArray<ndim, T> &tgt, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);

		/**
		 * A global array that is partitioned across a process group according to a domain
		 */
		template<int ndim,typename T>
		class GlobalArray
		{
			friend class Access<ndim,T>;
                        friend void get<>(GlobalArray<ndim, T> &tgt, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);
                        friend void put<>(GlobalArray<ndim, T> &tgt, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);
                        friend void accumulate<>(GlobalArray<ndim, T> &tgt, coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);
	
		public:
                        typedef T value_type;
			typedef GARef<ndim,T> storage;
			typedef Access<ndim,T> accessor;
			static const /*constexpr*/ int number_of_dimensions = ndim;
			typedef std::array<Boundary<ndim,T>,ndim> bounds;

		private:
			// before allocate
			const Domain<ndim>* dom;
			coords<ndim> gw;
			bool gc;
			bounds bd;

			// after allocate
		public:
			T* ptr;
			std::unique_ptr<GAImpl<ndim, T>> impl;

                private:
			coords<ndim+1> ld;

                        // ghost regions
			std::array<coords_range<ndim>,ndim> ghost_back;
			std::array<coords_range<ndim>,ndim> ghost_front;

			// extra
			mutable int lc;

			void allocate();
			void deallocate();
			void reset();

                public:
                        std::array<GABuf<ndim,T>, ndim> g_out_back, g_out_front, g_in_back, g_in_front;

			INLINE T& da(coords<ndim> i) const;
			INLINE T* pa(coords<ndim> i) const;
                        INLINE coord offset(coords<ndim> i) const;
                        INLINE coord byte_offset(coords<ndim> i) const;

			std::ostream& log_out() const;

			/**
			 * The domain
			 */
			INLINE const Domain<ndim>& domain() const;

			/**
			 * The number of ghost cells for each dimension.
			 */
			INLINE coords<ndim> ghost_width() const;

			/**
			 * Whether corner ghost cells are maintained.
			 */
			INLINE bool ghost_corners() const;

			/**
			 * Returns whether this global array is active.
			 */
			INLINE /*explicit*/ operator bool() const;

			/**
			 * The total/inner/outer region of elements
			 */
			INLINE coords_range<ndim> region() const;
			INLINE coords_range<ndim> local() const;
			INLINE coords_range<ndim> inner() const;
			INLINE coords_range<ndim> outer_front(int) const;
			INLINE coords_range<ndim> outer_back(int) const;
                        INLINE std::vector<coords_range<ndim>> outer() const;
                        INLINE std::vector<coords_range<ndim>> split(coords_range<ndim>) const;

			INLINE coords_range<ndim> ghost_in_front(int) const;
			INLINE coords_range<ndim> ghost_in_back(int) const;
			INLINE coords_range<ndim> ghost_out_front(int) const;
			INLINE coords_range<ndim> ghost_out_back(int) const;

			/**
			 * Select a region of elements as destination
			 */
			GADest<ndim,T> region(coords_range<ndim>, bool = true);

			/**
			 * Construct a GlobalArray (collective).
			 * The global array will not be active until it is assigned.
			 */
			GlobalArray();

			/**
			 * Construct a GlobalArray (collective).
			 * @param domain the domain for the distribution of data
			 * @param ghost_width number of ghost cells for each dimension
			 * @param ghost_corners whether to maintain the corner ghost cells
			 */
			GlobalArray(const Domain<ndim>& domain, coords<ndim> ghost_width = coords<ndim>(), bool ghost_corners = false, bounds bd = bounds());

			GlobalArray(const GlobalArray<ndim,T>& other, bool copy);

			/**
			 * Destruct a GlobalArray (collective).
			 * If active, the memory of the global array is released
			 */
			~GlobalArray();

			// Move semantics
			GlobalArray(const GlobalArray<ndim,T>& other) = delete;

			GlobalArray(GlobalArray<ndim,T>&& other);
			GlobalArray<ndim,T>& operator=(GlobalArray<ndim,T>&& other);

			// Copy from source
			GlobalArray<ndim,T>& operator=(const GlobalArray<ndim,T>& other);
			template<typename S>
			GlobalArray<ndim,T>& operator=(const S&);
			template<typename S>
			GlobalArray<ndim,T>& operator<<=(const S&);

			/**
			 * Update the ghost cells with data from their original locations (collective).
			 * RMA operations cannot overlap with local access.
			 *
			 * iupdate is a non-blocking variant if circumstances permit this (aync communication
			 * support, no ghost corners)
			 *
			 * @param k update phase (only used for general boundaries)
			 */
			void update(long k=0) const;
			Future<void> iupdate(long k=0) const;

			/**
			 * Get remote range (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to retrieve from
			 * @param ld the strides to use for buf (default: determined by range)
			 * @param buf the target buffer
			 */
                        void get(coords_range<ndim> range_src, GlobalArray<ndim, T> &src, coords_range<ndim> range_tgt);
			void get(coords_range<ndim> range, T* buf) const;
			void get(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, T* buf) const;

			/**
			 * Put remote range (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to store to
			 * @param ld the strides to use for buf (default: determined by range)
			 * @param buf the source buffer
			 */
			void put(coords_range<ndim> range_tgt, GlobalArray<ndim, T> &src, coords_range<ndim> range_src);
			void put(coords_range<ndim> range, const T* buf);
			void put(coords_range<ndim> range, std::array<std::size_t,ndim-1> ld, const T* buf);

			/**
			 * Increase remote range with local values (one-sided).
			 * RMA operations cannot overlap with local access.
			 * @param range the area of the global array to update
			 * @param ld the strides to use for buf (default: determined by range)
			 * @param buf the source buffer
			 */
			void accumulate(coords_range<ndim> range_tgt, GlobalArray<ndim, T> &src, coords_range<ndim> range_src);

			/**
			 * Dump values of global array to a file
			 */
			template<typename = void>
			void dump(coords_range<ndim> range, std::string filename);
		};

                template<int ndim,typename T>
                std::ostream& operator<<(std::ostream& out, const GlobalArray<ndim,T> &ga);

		template<int ndim, typename T>
		class GADest {
			friend class GlobalArray<ndim,T>;
			GlobalArray<ndim,T>& ga;
                        coords_range<ndim> region;
                        bool pack;
			GADest(GlobalArray<ndim,T>& ga, coords_range<ndim> r, bool pack);
		public:
			~GADest();
			GADest(const GADest<ndim,T>& gad);
			GADest& operator=(const GADest<ndim,T>& gad) = delete;
			template<typename S>
			GADest<ndim,T>& operator=(const S& src);
		};

		template<int ndim, typename T>
		class GARef {
			const GlobalArray<ndim,T>& ga;
		public:
			GARef(const GlobalArray<ndim,T>& ga);
			~GARef();
			INLINE operator const GlobalArray<ndim,T>&() const;
			INLINE const Domain<ndim>& domain() const;
			INLINE coords_range<ndim> region() const;
		};

		template<int ndim, typename T>
		struct GABuf {
                        coords_range<ndim> r;
                        T* ptr;
                        mutable bool full;

			GABuf() : full(false) {}
			GABuf(coords_range<ndim> r, T* ptr) : r(r), ptr(ptr), full(false) {}

                        T& da(coords<ndim> i) const {
                                auto off = (i - r.lower).offset(r.stride());
                                assert(off < r.count());
                                return ptr[off];
                        }
		};

		// Inline function implementations

		template<int ndim, typename T>
		inline const Domain<ndim>& GlobalArray<ndim,T>::domain() const {
			return *dom;
		}

		template<int ndim, typename T>
		inline coords<ndim> GlobalArray<ndim,T>::ghost_width() const {
			return gw;
		}

		template<int ndim, typename T>
		inline bool GlobalArray<ndim,T>::ghost_corners() const {
			return gc;
		}

		template<int ndim, typename T>
		inline GlobalArray<ndim,T>::operator bool() const {
			return dom != 0;
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::region() const {
			return domain().total();
		}

		template<int ndim, typename T>
		inline T* GlobalArray<ndim,T>::pa(coords<ndim> i) const {
			return ptr + offset(i);
                }

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::local() const {
			return domain().local();
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::inner() const {
			auto r = domain().local();
                        r.lower += ghost_width();
                        r.upper -= ghost_width();
                        return r;
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::outer_front(int d) const {
                        assert(d >= 0 && d < ndim);
                        const auto gw = ghost_width();
			auto r = domain().local();
                        for(int i = 0; i<d; ++i) { 
                            r.lower[i] += gw[i];
                            r.upper[i] -= gw[i];
                        }
                        r.upper[d] = r.lower[d] + ghost_width()[d];
                        return r;
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::outer_back(int d) const {
                        assert(d >= 0 && d < ndim);
                        const auto gw = ghost_width();
                        auto r = domain().local();
                        for(int i = 0; i<d; ++i) {
                            r.lower[i] += gw[i];
                            r.upper[i] -= gw[i];
                        }
                        r.lower[d] = r.upper[d] - gw[d];
                        return r;
		}

		template<int ndim, typename T>
		inline std::vector<coords_range<ndim>> GlobalArray<ndim,T>::outer() const {
                        auto r = std::vector<coords_range<ndim>>(2*ndim);
                        seq<0,ndim>::for_each([this,&r](int d) { 
                                r[2*d  ] = outer_front(d);
                                r[2*d+1] = outer_back(d);
                        });
                        return r;
		}

                template<int ndim, typename T>
                inline std::vector<coords_range<ndim>> GlobalArray<ndim,T>::split(coords_range<ndim> r) const {
                    int count = 0;
                    if (r == region()) {
                        //log_out() << "split " << r << " into " << std::endl;
                        auto ret = outer();
                        ret.push_back(inner());
                        for(auto r : ret) {
                            //log_out() << " " << r << std::endl;
                            count += r.count();
                        }
                        assert(count == local().count());
                        return ret;
                    }
                    return std::vector<coords_range<ndim>>(1,r);
                }


		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::ghost_in_front(int d) const {
                    return ghost_front[d];
                }

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::ghost_in_back(int d) const {
                    return ghost_back[d];
                }

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::ghost_out_front(int d) const {
                    return ghost_front[d].adj(d, -1);
                }

		template<int ndim, typename T>
		inline coords_range<ndim> GlobalArray<ndim,T>::ghost_out_back(int d) const {
                    return ghost_back[d].adj(d, 1);
                }

		template<int ndim, typename T>
		inline T& GlobalArray<ndim,T>::da(coords<ndim> i) const {
			return ptr[offset(i)];
		}

       		template<int ndim, typename T>
		inline coord GlobalArray<ndim,T>::offset(coords<ndim> i) const {
			return (i + gw).offset(ld);
		}

       		template<int ndim, typename T>
		inline coord GlobalArray<ndim,T>::byte_offset(coords<ndim> i) const {
			return offset(i) * sizeof(T);
		}

		template<int ndim, typename T>
		inline GARef<ndim,T>::operator const GlobalArray<ndim,T>&() const {
			return ga;
		}

		template<int ndim, typename T>
		inline const Domain<ndim>& GARef<ndim,T>::domain() const {
			return ga.domain();
		}

		template<int ndim, typename T>
		inline coords_range<ndim> GARef<ndim,T>::region() const {
			return ga.region();
		}

		// Generic members

		template<int ndim, typename T> template<typename S>
		GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator=(const S& src) {
			region(domain().total(), false) = src;
			return *this;
		}

		template<int ndim, typename T> template<typename S>
		GlobalArray<ndim,T>& GlobalArray<ndim,T>::operator<<=(const S& src) {
			region(domain().total(), true) = src;
			return *this;
		}

		template<int ndim, typename T> template<typename S>
		GADest<ndim,T>& GADest<ndim,T>::operator=(const S& src) {
			static_assert(S::number_of_dimensions == ndim, "source dimensionality");
			assert(ga.domain() == src.domain());
			Access<ndim,T> d(ga);
			const typename S::accessor s(src);

                        if (pack) {
                            auto regions = ga.split(region);
                            for(auto r : regions) {
                                assert(src.region().contains(r));
                                std::vector<GABuf<ndim,T> *> ghosts;
                                GABuf<ndim,T> *full_ghost = 0;

                                auto push = [this, &full_ghost, &ghosts,&r](GABuf<ndim,T> *g) {
                                    auto overlap = g->r.overlap(r); 
                                    if (g->r.count() > 0) ga.log_out() << " overlap " << 100 * overlap.count() / g->r.count() << "%" << std::endl;
                                    if (!full_ghost && overlap.count() > 0 && overlap.count() == g->r.count()) {
                                        full_ghost = g;
                                        g->full = true;
                                    } else if (overlap.count() > 0) {
                                        ghosts.push_back(g);
                                        g->full = true;
                                    }
                                };

                                seq<0,ndim>::for_each([this, &r, &full_ghost, &ghosts, &push](int d) {
                                        push(&ga.g_out_front[d]);
                                        push(&ga.g_out_back[d]);
                                        });

                                //ga.log_out() << "total ghosts that overlap: " << ghosts.size() << " out of " <<  regions.size() << std::endl;

                                if(full_ghost && ghosts.size()) {
                                    //SHARK_COUNTER("pack_full+ghost");
                                    ga.domain().for_each(r, [this, &d, &s, &full_ghost, &ghosts](coords<ndim> i){
                                            auto v = d(i) = s(i);
                                            full_ghost->da(i) = v;
                                            for(auto g : ghosts) if (g->r.contains(i)) g->da(i) = v;
                                    });
                                } else if (full_ghost) {
                                    //SHARK_COUNTER("pack_onlyfull");
                                    ga.domain().for_each(r, [this, &d, &s, &full_ghost](coords<ndim> i){
                                            auto v = d(i) = s(i);
                                            full_ghost->da(i) = v;
                                    });
                                } else if (ghosts.size()) {
                                    //SHARK_COUNTER("pack_ghost");
                                    ga.domain().for_each(r, [this, &d, &s, &ghosts](coords<ndim> i){
                                            auto v = d(i) = s(i);
                                            for(auto g : ghosts) if (g->r.contains(i)) g->da(i) = v;
                                    });
                                } else {
                                    //SHARK_COUNTER("pack_no_ghost");
                                    ga.domain().for_each(r, [this, &d, &s](coords<ndim> i){
                                            d(i) = s(i);
                                    });
                                }
                            }
                        } else  {
                            //SHARK_COUNTER("pack_no_pack");
                            ga.domain().for_each(region, [this, &d, &s](coords<ndim> i){
                                    d(i) = s(i);
                            });
                        }

                        return *this;
                }

	}

}

#endif
