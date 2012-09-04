#ifndef __SHARK_COMMON_HPP
#define __SHARK_COMMON_HPP

#include <type_traits>                 // std::enable_if

#define INLINE inline __attribute((__no_instrument_function__))

namespace shark {

	/**
	 * Type of a single coordinate
	 */
	typedef long coord;

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
		class GlobalArray;
		template<int,typename>
		class Access;

	}

	struct test_result {
		long fails;
		long checks;
	};

}

#endif
