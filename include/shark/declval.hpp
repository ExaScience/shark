/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#ifndef __SHARK_DECLVAL_HPP
#define __SHARK_DECLVAL_HPP

#include <utility>        // std::declval

#ifdef SHARK_GLIBCXX_TR1IMPL
// Workaround for std::declval
#include <type_traits>
namespace std {
	template<typename T>
	typename add_rvalue_reference<T>::type declval();
}
#endif

#endif
