#ifndef __SHARK_ABS_HPP
#define __SHARK_ABS_HPP

#include <cstdlib>  // std::abs(int)
#include <cmath>    // std::abs(float)

#include "common.hpp"

namespace shark {
	
	namespace ndim {

		// abs function for builtin types like int, float, double,... where argument-dependent lookup does not work
		using std::abs;

	}

}

#endif
