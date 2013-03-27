/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_COMMON_HPP
#define __SHARK_COMMON_HPP

#include <type_traits>                 // std::enable_if
#include <bitset>                      // std::bitset
#include <ostream>                     // std::ostream
#if defined(SHARK_PTHREAD_SCHED)
#include <functional>                  // std::function
#endif

/*
 * Experimental configuration options not managed by configure
 */
//#define SHARK_RANGE
//#define SHARK_THREAD_BLOCK_DIST
//#define SHARK_OMP_TDIST

#define INLINE inline __attribute((__no_instrument_function__))

namespace shark {

	/**
	 * Type of a single coordinate
	 */
	typedef long coord;

	static const int verbose_alloc = 0;
	static const int verbose_update = 1;
	static const int verbose_rma = 2;
	static const int verbose_collective = 3;
	static const int verbose_end = 4;
	typedef std::bitset<verbose_end> verbosity_mask;

	/**
	 * Log mask to configure verbosity of logging of operations.
	 */
	extern verbosity_mask log_mask;

	/**
	 * Pointer to the stream used for log output.
	 * Should be non-NULL when the log_mask indicates output.
	 */
	extern std::ostream* log_out;

#if defined(SHARK_PTHREAD_SCHED)
	void ThreadWork(std::function<void(int)> w);
#endif

	// (Available only to shark)
	namespace {

		// Empty function to signify unused named parameters
		template<typename... Ts>
		INLINE void unused(const Ts&...);

		template<typename... Ts>
		inline void unused(const Ts&...) { }

		/**
		 * Compile-time sequence consisting of [start, start+1, ..., start+count)
		 */
		template<int start, int count>
		class seq {
		public:
			template<typename Func>
			INLINE static void for_each(const Func& f);
			template<typename Func>
			INLINE static bool all_of(const Func& f);
			template<typename Func, typename T>
			INLINE static T product(const Func& f, const T& base);
			template<typename Func, typename T>
			INLINE static T sum(const Func& f, const T& base);
		};
		template<int start>
		class seq<start,0> {
		public:
			template<typename Func>
			INLINE static void for_each(const Func& f);
			template<typename Func>
			INLINE static bool all_of(const Func& f);
			template<typename Func, typename T>
			INLINE static T product(const Func& f, const T& base);
			template<typename Func, typename T>
			INLINE static T sum(const Func& f, const T& base);
		};


		// Inline function implementations
		
		template<int start, int count> template<typename Func>
		inline void seq<start,count>::for_each(const Func& f) {
			f(start);
			seq<start+1,count-1>::for_each(f);
		}

		template<int start> template<typename Func>
		inline void seq<start,0>::for_each(const Func&) {
			// Done
		}

		template<int start, int count> template<typename Func>
		inline bool seq<start,count>::all_of(const Func& f) {
			return f(start) && seq<start+1,count-1>::all_of(f);
		}

		template<int start> template<typename Func>
		inline bool seq<start,0>::all_of(const Func&) {
			// Base case
			return true;
		}

		template<int start, int count> template<typename Func, typename T>
		inline T seq<start,count>::product(const Func& f, const T& base) {
			return f(start) * seq<start+1,count-1>::product(f, base);
		}

		template<int start> template<typename Func, typename T>
		inline T seq<start,0>::product(const Func&, const T& base) {
			// Base case
			return base;
		}

		template<int start, int count> template<typename Func, typename T>
		inline T seq<start,count>::sum(const Func& f, const T& base) {
			return f(start) + seq<start+1,count-1>::sum(f, base);
		}

		template<int start> template<typename Func, typename T>
		inline T seq<start,0>::sum(const Func&, const T& base) {
			// Base case
			return base;
		}
	}

	// Forward declaration
	class Group;

	namespace ndim {

		// Forward declarations
		template<int>
		struct coords;
		template<int>
		struct coords_range;
		template<int ndim>
		class Domain;
		template<int,typename>
		class SparseArray;
		template<int,typename>
		class Boundary;
		template<int,typename>
		class GlobalArray;
		template<int,typename>
		class Access;

	}

	struct test_result {
		long fails;
		long checks;

		INLINE test_result& operator+=(const test_result&);
	};

	inline test_result& test_result::operator+=(const test_result& other) {
		fails += other.fails;
		checks += other.checks;
		return *this;
	}

}

#endif
