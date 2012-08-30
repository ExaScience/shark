
#include <shark/domain.hpp>
#include <shark/sparsearray.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

template<int ndim,typename T>
unsigned short SparseArray<ndim,T>::full_cw(const Domain<ndim>& dom) {
	unsigned short c = 0;
	for(int d = 0; d < ndim; d++)
		while(dom.n[d] >> c != 0)
			++c;
	return c;
}

template<int ndim,typename T>
coords<ndim+1> SparseArray<ndim,T>::stride(const Domain<ndim>& dom) {
	coords<ndim+1> ld;
	ld[ndim] = 1;
	for(int d = ndim-1; d >= 0; d--)
		ld[d] = ld[d+1] * dom.n[d];
	return ld;
}

template<int ndim,typename T>
coords<ndim+1> SparseArray<ndim,T>::coarse_stride(const Domain<ndim>& dom, unsigned short cw) {
	coords<ndim+1> ld;
	ld[ndim] = 1;
	for(int d = ndim-1; d >= 0; d--)
		ld[d] = ld[d+1] * (((dom.n[d] - 1) >> cw) + 1);
	return ld;
}

template<int ndim,typename T>
SparseArray<ndim,T>::SparseArray(const Domain<ndim>& dom, unsigned short cw):
	dom(dom), cw(cw), ld(stride(dom)), cld(coarse_stride(dom,cw)), ptr(new T[ld[0]]), coarse(cld[0]) {
	init();
}

template<int ndim,typename T>
SparseArray<ndim,T>::SparseArray(const Domain<ndim>& dom):
	dom(dom), cw(full_cw(dom)), ld(stride(dom)), cld(coarse_stride(dom,cw)), ptr(new T[ld[0]]), coarse(cld[0]) {
	init();
}

template<int ndim,typename T>
SparseArray<ndim,T>::SparseArray(const SparseArray<ndim,T>& other):
	dom(other.dom), cw(other.cw), ld(other.ld), cld(other.cld), ptr(new T[ld[0]]), coarse(other.coarse) {
	dom.total().for_each([this,&other](coords<ndim> i) {
		coord off = i.offset(ld);
		ptr[off] = other.ptr[off];
	});
}

template<int ndim,typename T>
SparseArray<ndim,T>& SparseArray<ndim,T>::operator=(const SparseArray<ndim,T>& other) {
	assert(dom == other.dom);
	assert(cw == other.cw);
	coarse = other.coarse;
	dom.total().for_each([this,&other](coords<ndim> i) {
		coord off = i.offset(ld);
		ptr[off] = other.ptr[off];
	});
	return *this;
}

template<int ndim,typename T>
SparseArray<ndim,T>::~SparseArray() {
	delete[] ptr;
}

template<int ndim,typename T>
void SparseArray<ndim,T>::init() {
	dom.total().for_each([this](coords<ndim> i) {
		ptr[i.offset(ld)] = T();
	});
}

template<int ndim,typename T>
void SparseArray<ndim,T>::clear() {
	init();
	coarse.assign(cld[0], false);
}

template<int ndim,typename T>
void SparseArray<ndim,T>::operator+=(const SparseArray<ndim,T>& other) {
	assert(dom == other.dom);
	assert(cw == other.cw);

	for(coord k = 0; k < cld[0]; k++)
		if(other.coarse[k]) {
			coarse[k] = true;

			coords_range<ndim> r;
			coord kk = k;
			for(int d = 0; d < ndim; d++) {
				coord fk = kk / cld[d+1];
				r.lower[d] = fk << cw;
				r.upper[d] = min((fk+1) << cw, dom.n[d]);
				kk = kk % cld[d+1];
			}
			r.for_each([this,&other](coords<ndim> i) {
				coord off = i.offset(ld);
				ptr[off] += other.ptr[off];
			});
		}
}


// Set-up instantiations

#include "types"

#define SYMBDT(d,T) template class shark::ndim::SparseArray<d,T>; 
#include "inst_dimtype"
#undef SYMBDT

