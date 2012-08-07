
#include <shark/coords_range.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim>
ostream& shark::ndim::operator<<(ostream& out, const coords_range<ndim>& r) {
	out << "[" << r.lower[0] << ":" << r.upper[0] << ")";
	for_each<1,ndim>([&out,&r](int d) { out << " x " << "[" << r.lower[d] << ":" << r.upper[d] << ")"; });
	return out;
}


// Set-up instantiations

#define SYMB(d) template ostream& shark::ndim::operator<< <d>(ostream&, const coords_range<d>&);
#include "inst_dim"
#undef SYMB
