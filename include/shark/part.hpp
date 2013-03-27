/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_PART_HPP
#define __SHARK_PART_HPP

#include "vec.hpp"

namespace shark {
	
	namespace ndim {

		template<int ndim, typename T>
		struct part {
			vec<ndim,T> x;
			vec<ndim,T> v;
			T w;
		};

	}

}

#endif
