
#include <shark/region.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim>
Region<ndim>::Region(const Domain<ndim>& domain)
	: domain(domain)
	, r(domain.local())
{ }

template<int ndim>
Region<ndim>::Region(const Domain<ndim>& domain, coords_range<ndim> r)
	: domain(domain)
	// Use comma operator to invoke overlap
	, r((r = domain.local().overlap(r),r))
{ }

template<int ndim>
Region<ndim>::~Region() { }

template<int ndim>
bool Region<ndim>::hasData() const {
	return seq<0,ndim>::all_of([this](int d) { return r.lower[d] < r.upper[d]; } );
}

// Set-up instantiations

#define SYMB(d) template class Region<d>;
#include "inst_dim"
#undef SYMB

