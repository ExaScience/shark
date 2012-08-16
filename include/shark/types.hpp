#ifndef __SHARK_TYPES_HPP
#define __SHARK_TYPES_HPP

#include "common.hpp"

namespace shark {

	namespace types2d {
		typedef shark::ndim::coords<2> coords;
		typedef shark::ndim::coords_range<2> coords_range;
		typedef shark::ndim::Domain<2> Domain;
		typedef shark::ndim::GlobalArray<2,int> GlobalArrayI;
		typedef shark::ndim::GlobalArray<2,long> GlobalArrayL;
		typedef shark::ndim::GlobalArray<2,float> GlobalArrayF;
		typedef shark::ndim::GlobalArray<2,double> GlobalArrayD;
		typedef shark::ndim::Access<2,int> AccessI;
		typedef shark::ndim::Access<2,long> AccessL;
		typedef shark::ndim::Access<2,float> AccessF;
		typedef shark::ndim::Access<2,double> AccessD;
		typedef shark::ndim::vec<2,double> vecD;
	}

	namespace types3d {
		typedef shark::ndim::coords<3> coords;
		typedef shark::ndim::coords_range<3> coords_range;
		typedef shark::ndim::Domain<3> Domain;
		typedef shark::ndim::GlobalArray<3,int> GlobalArrayI;
		typedef shark::ndim::GlobalArray<3,long> GlobalArrayL;
		typedef shark::ndim::GlobalArray<3,float> GlobalArrayF;
		typedef shark::ndim::GlobalArray<3,double> GlobalArrayD;
		typedef shark::ndim::Access<3,int> AccessI;
		typedef shark::ndim::Access<3,long> AccessL;
		typedef shark::ndim::Access<3,float> AccessF;
		typedef shark::ndim::Access<3,double> AccessD;
		typedef shark::ndim::vec<3,double> vecD;
	}

}

#endif
