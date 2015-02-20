/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_PART_HPP
#define __SHARK_PART_HPP

#include "vec.hpp"

namespace shark {

	namespace ndim {

		struct part_position {
			int pos;
			float off;
		};

		template<int ndim>
		struct part {
			vec<ndim,part_position> x;
			vec<ndim,double> v;
			double w;
		};

	}

}

#endif
