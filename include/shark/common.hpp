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

		/**
		 * Template functions to iterate over a compile-time range [begin..end)
		 */
		template<int begin, int end, typename Func>
		INLINE typename std::enable_if<begin < end>::type for_each(const Func& f);
		template<int begin, int end, typename Func>
		INLINE typename std::enable_if<begin == end>::type for_each(const Func& f);
		template<int begin, int end, typename Func>
		INLINE typename std::enable_if<begin < end,bool>::type all_of(const Func& f);
		template<int begin, int end, typename Func>
		INLINE typename std::enable_if<begin == end,bool>::type all_of(const Func& f);
		template<int begin, int end, typename Func, typename T>
		INLINE typename std::enable_if<begin < end,T>::type product(const Func& f, const T& identity);
		template<int begin, int end, typename Func, typename T>
		INLINE typename std::enable_if<begin == end,T>::type product(const Func& f, const T& identity);
		template<int begin, int end, typename Func, typename T>
		INLINE typename std::enable_if<begin < end,T>::type sum(const Func& f, const T& zero);
		template<int begin, int end, typename Func, typename T>
		INLINE typename std::enable_if<begin == end,T>::type sum(const Func& f, const T& zero);

		// Inline function implementations

		template<int begin, int end, typename Func>
		inline typename std::enable_if<begin < end>::type for_each(const Func& f) {
			f(begin);
			for_each<begin+1,end,Func>(f);
		}

		template<int begin, int end, typename Func>
		inline typename std::enable_if<begin == end>::type for_each(const Func&) {
			// Done
		}

		template<int begin, int end, typename Func>
		inline typename std::enable_if<begin < end,bool>::type all_of(const Func& f) {
			return f(begin) && all_of<begin+1,end,Func>(f);
		}

		template<int begin, int end, typename Func>
		inline typename std::enable_if<begin == end,bool>::type all_of(const Func&) {
			// Base case
			return true;
		}

		template<int begin, int end, typename Func, typename T>
		inline typename std::enable_if<begin < end,T>::type product(const Func& f, const T& identity) {
			return f(begin) * product<begin+1,end,Func,T>(f, identity);
		}

		template<int begin, int end, typename Func, typename T>
		inline typename std::enable_if<begin == end,T>::type product(const Func&, const T& identity) {
			// Base case
			return identity;
		}

		template<int begin, int end, typename Func, typename T>
		inline typename std::enable_if<begin < end,T>::type sum(const Func& f, const T& zero) {
			return f(begin) + sum<begin+1,end,Func,T>(f, zero);
		}

		template<int begin, int end, typename Func, typename T>
		inline typename std::enable_if<begin == end,T>::type sum(const Func&, const T& zero) {
			// Base case
			return zero;
		}
	}

	// Forward declaration
	class Group;

	namespace ndim {

		// Forward declarations
		template<int ndim>
		class Domain;
		template<int>
		class Region;
		template<int,typename>
		class SparseArray;
		template<int,typename>
		class GlobalArray;
		template<int,typename>
		class Access;
	}

}

#endif
