
#include <shark/coords.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim>
ostream& shark::ndim::operator<<(ostream& out, const coords<ndim>& i) {
	out << i[0];
	seq<1,ndim-1>::for_each([&out,&i](int d) { out << "," << i[d]; });
	return out;
}

// Set-up instantiations

#define SYMB(d) template ostream& shark::ndim::operator<< <d>(ostream&, const coords<d>&);
#include "inst_dim"
#undef SYMB

