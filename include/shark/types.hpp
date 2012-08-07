#ifndef __SHARK_TYPES_HPP
#define __SHARK_TYPES_HPP

#include "common.hpp"

namespace shark {

	namespace types3d {
		typedef shark::ndim::coords<3> coords;
		typedef shark::ndim::coords_range<3> coords_range;
		typedef shark::ndim::Domain<3> Domain;
		typedef shark::ndim::Region<3> Region;
		typedef shark::ndim::GlobalArray<3,int> GlobalArrayI;
		typedef shark::ndim::GlobalArray<3,long> GlobalArrayL;
		typedef shark::ndim::GlobalArray<3,float> GlobalArrayF;
		typedef shark::ndim::GlobalArray<3,double> GlobalArrayD;
		typedef shark::ndim::Access<3,int> AccessI;
		typedef shark::ndim::Access<3,long> AccessL;
		typedef shark::ndim::Access<3,float> AccessF;
		typedef shark::ndim::Access<3,double> AccessD;
	}

}

#endif
