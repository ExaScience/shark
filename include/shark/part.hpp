/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
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

                        part_position& operator+=(const part_position& other) {
                            pos += other.pos;
                            off += other.off;
                            return *this;
                        }

                        part_position operator+(const part_position& other) const {
                            part_position r;
                            r.pos = pos + other.pos;
                            r.off = off + other.off;
                            return r;
                        }
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

                        inline part<ndim>& operator+=(const part<ndim>& other) {
                                x += other.x;
                                v += other.v;
                                w += other.w;
                                return *this;
                        }

                        inline part<ndim> operator+(const part<ndim>& other) const {
                                part<ndim> r;
                                r.x = x + other.x;
                                r.v = v + other.v;
                                r.w = w + other.w;
                                return r;
                        }
		};

                template<int ndim>
       		inline std::ostream & operator<<(std::ostream &os, const part<ndim> &v) {
                        os << "[" << v.x << ", " << v.v << ", " << v.w << "]";;
			return os;
		}

                
 	}

}

#endif
