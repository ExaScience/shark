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

	  inline bool operator<= (const part_position& p1, const part_position& p2) {
	    if (p1.pos < p2.pos) return true;
	    else if (p1.pos == p2.pos) return p1.off <= p2.off;
	    else return false;
	  }

		template<int ndim>
		struct part {
			vec<ndim,part_position> x;
			vec<ndim,double> v;
			double w;
		};

	}

}

#endif
