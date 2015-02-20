/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __SHARK_PART_HPP
#define __SHARK_PART_HPP

#include "vec.hpp"

namespace shark {

	namespace ndim {

		/**
		 * Position of a particle as a cell number and offset
		 */
		struct part_position {
			int pos;
			float off;
		};

		INLINE bool operator<=(const part_position& p1, const part_position& p2);

		// Inline members

		inline bool operator<=(const part_position& p1, const part_position& p2) {
			return p1.pos < p2.pos || (p1.pos == p2.pos && p1.off <= p2.off);
		}

		/**
		 * Particles with position, velocity, charge
		 */
		template<int ndim>
		struct part {
			vec<ndim,part_position> x;
			vec<ndim,double> v;
			double w;
		};

	}

}

#endif
