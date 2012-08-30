
#include <shark/coords_range.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim>
coords_range<ndim> coords_range<ndim>::overlap(coords_range<ndim> other) const {
	coords_range<ndim> r;
	seq<0,ndim>::for_each([this, &r, &other](int d) {
		r.lower[d] = std::max(lower[d], other.lower[d]);
		r.upper[d] = std::min(upper[d], other.upper[d]);
	});
	return r;
}

template<int ndim>
bool coords_range<ndim>::contains(coords<ndim> i) const {
	return seq<0,ndim>::all_of([this, &i](int d) {
		return lower[d] <= i[d] && upper[d] > i[d];
	});
}

template<int ndim>
bool coords_range<ndim>::contains(coords_range<ndim> other) const {
	return seq<0,ndim>::all_of([this, &other](int d) {
		return lower[d] <= other.lower[d] && upper[d] >= other.upper[d];
	});
}

template<int ndim>
ostream& shark::ndim::operator<<(ostream& out, const coords_range<ndim>& r) {
	out << "[" << r.lower[0] << ":" << r.upper[0] << ")";
	seq<1,ndim-1>::for_each([&out,&r](int d) { out << " x " << "[" << r.lower[d] << ":" << r.upper[d] << ")"; });
	return out;
}

// Set-up instantiations

#define SYMBD(d) template struct shark::ndim::coords_range<d>;
#include "inst_dim"
#undef SYMBD

#define SYMBD(d) template ostream& shark::ndim::operator<< <d>(ostream&, const coords_range<d>&);
#include "inst_dim"
#undef SYMBD
