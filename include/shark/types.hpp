#ifndef __SHARK_TYPES_HPP
#define __SHARK_TYPES_HPP

#include "common.hpp"

namespace shark {

	namespace types2d {
		typedef shark::ndim::coords<2> coords;
		typedef shark::ndim::coords_range<2> coords_range;
		typedef shark::ndim::Domain<2> Domain;
		typedef shark::ndim::vec<2,double> vecD;
		typedef shark::ndim::part<2,double> partD;
		typedef shark::ndim::GlobalArray<2,int> GlobalArrayI;
		typedef shark::ndim::GlobalArray<2,long> GlobalArrayL;
		typedef shark::ndim::GlobalArray<2,float> GlobalArrayF;
		typedef shark::ndim::GlobalArray<2,double> GlobalArrayD;
		typedef shark::ndim::GlobalArray<2,vecD> GlobalArrayV;
		typedef shark::ndim::GlobalArray<2,partD> GlobalArrayP;
		typedef shark::ndim::Access<2,int> AccessI;
		typedef shark::ndim::Access<2,long> AccessL;
		typedef shark::ndim::Access<2,float> AccessF;
		typedef shark::ndim::Access<2,double> AccessD;
		typedef shark::ndim::Access<2,vecD> AccessV;
		typedef shark::ndim::Access<2,partD> AccessP;
		typedef shark::ndim::Boundary<2,int> BoundaryI;
		typedef shark::ndim::Boundary<2,long> BoundaryL;
		typedef shark::ndim::Boundary<2,float> BoundaryF;
		typedef shark::ndim::Boundary<2,double> BoundaryD;
		typedef shark::ndim::Boundary<2,vecD> BoundaryV;
		typedef shark::ndim::Boundary<2,partD> BoundaryP;
		typedef shark::ndim::SparseArray<2,int> SparseArrayI;
		typedef shark::ndim::SparseArray<2,long> SparseArrayL;
		typedef shark::ndim::SparseArray<2,float> SparseArrayF;
		typedef shark::ndim::SparseArray<2,double> SparseArrayD;
		typedef shark::ndim::SparseArray<2,vecD> SparseArrayV;

		// Hack: Allow ADL for function call with explicit template argument according to 14.8.1/8
		template<typename T> void coord_val();
	}

	namespace types3d {
		typedef shark::ndim::coords<3> coords;
		typedef shark::ndim::coords_range<3> coords_range;
		typedef shark::ndim::Domain<3> Domain;
		typedef shark::ndim::vec<3,double> vecD;
		typedef shark::ndim::part<3,double> partD;
		typedef shark::ndim::GlobalArray<3,int> GlobalArrayI;
		typedef shark::ndim::GlobalArray<3,long> GlobalArrayL;
		typedef shark::ndim::GlobalArray<3,float> GlobalArrayF;
		typedef shark::ndim::GlobalArray<3,double> GlobalArrayD;
		typedef shark::ndim::GlobalArray<3,vecD> GlobalArrayV;
		typedef shark::ndim::GlobalArray<3,partD> GlobalArrayP;
		typedef shark::ndim::Access<3,int> AccessI;
		typedef shark::ndim::Access<3,long> AccessL;
		typedef shark::ndim::Access<3,float> AccessF;
		typedef shark::ndim::Access<3,double> AccessD;
		typedef shark::ndim::Access<3,vecD> AccessV;
		typedef shark::ndim::Access<3,partD> AccessP;
		typedef shark::ndim::Boundary<3,int> BoundaryI;
		typedef shark::ndim::Boundary<3,long> BoundaryL;
		typedef shark::ndim::Boundary<3,float> BoundaryF;
		typedef shark::ndim::Boundary<3,double> BoundaryD;
		typedef shark::ndim::Boundary<3,vecD> BoundaryV;
		typedef shark::ndim::Boundary<3,partD> BoundaryP;
		typedef shark::ndim::SparseArray<3,int> SparseArrayI;
		typedef shark::ndim::SparseArray<3,long> SparseArrayL;
		typedef shark::ndim::SparseArray<3,float> SparseArrayF;
		typedef shark::ndim::SparseArray<3,double> SparseArrayD;
		typedef shark::ndim::SparseArray<3,vecD> SparseArrayV;

		// Hack: Allow ADL for function call with explicit template argument according to 14.8.1/8
		template<typename T> void coord_val();
	}

}

#endif
