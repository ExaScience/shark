/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_ABS_HPP
#define __SHARK_ABS_HPP

#include <cstdlib>      // std::abs(int)
#include <cmath>        // std::abs(float)

namespace shark {
	
	// ADL for builtin types looks in the top namespace
	// 
	// Unfortunately ::abs(int) and ::fabs(float) are under different names
	//
	// Hence, we must use std::abs, which does have std::abs(int) and std::abs(float)
	using std::abs;

	template<typename T>
	INLINE auto adl_abs(T t) -> decltype(abs(t));

	template<typename T>
	INLINE auto adl_sin(T t) -> decltype(sin(t));

	template<typename T>
	INLINE auto adl_cos(T t) -> decltype(cos(t));

	template<typename T>
	inline auto adl_abs(T t) -> decltype(abs(t)) {
		return abs(t);
	}

	template<typename T>
	inline auto adl_sin(T t) -> decltype(sin(t)) {
		return sin(t);
	}

	template<typename T>
	inline auto adl_cos(T t) -> decltype(cos(t)) {
		return cos(t);
	}

}

#endif
