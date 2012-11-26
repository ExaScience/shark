/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_DECLVAL_HPP
#define __SHARK_DECLVAL_HPP

#include <utility>        // std::declval

// Temporary fix to enable declval for libstdc++ 4.4
#if defined(__GLIBCXX__) && __GLIBCXX__ < 20100414
#include <type_traits>
namespace std {
	template<typename T>
	typename add_rvalue_reference<T>::type declval();
}
#endif

#endif
