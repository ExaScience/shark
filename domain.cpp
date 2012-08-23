
#include <cassert>                     // assert

#include <shark/domain.hpp>
#include "comm_impl.hpp"

using namespace std;
using namespace shark;
using namespace shark::ndim;

namespace {

	vector<int> factor(int n) {
		vector<int> divs;
		int k = 2;
		while(n > 1) {
			while(n % k > 0)
				k++;
			divs.push_back(k);
			n /= k;
		}
		return divs;
	}

	// This function is a reimplementation of MPI_Dims_create

	void assign_dims(int n, int ndim, int dims[]) {
		vector<bool> skip(ndim);
		for(int d = 0; d < ndim; d++)
			if(dims[d] != 0) {
				skip[d] = true;
				n /= dims[d];
			} else
				dims[d] = 1;
		vector<int> facs = factor(n);
		int d = -1;
		for(auto it = facs.crbegin(); it != facs.crend(); ++it) {
			do {
				d++;
				if(d == ndim)
					d = 0;
			} while(skip[d]);
			dims[d] *= *it;
		}
	}

}

template<int ndim>
void Domain<ndim>::adjustProcs(int nprocs, pcoords& np) {
	assign_dims(nprocs, ndim, np.data());
}

template<int ndim>
array<int,ndim+1> Domain<ndim>::base(pcoords np) {
	array<int,ndim+1> b;
	b[ndim] = 1;
	for(int d = ndim-1; d >= 0; d--)
		b[d] = np[d] * b[d+1];
	return b;
}

template<int ndim>
typename Domain<ndim>::dists Domain<ndim>::distribution(coords<ndim> n, pcoords np) {
	dists nd = {{}};

	seq<0,ndim>::for_each([&nd, &n, &np](int d) {
		auto& dist(nd[d]);
		dist.resize(np[d]+1);

		int div = n[d] / np[d];
		int rem = n[d] % np[d];
		dist[0] = 0;
		for(int k = 0; k < np[d]; k++)
			dist[k+1] = dist[k] + div + (k < rem ? 1 : 0);

	});

	return nd;
}

template<int ndim>
bool Domain<ndim>::consistentDistribution() const {
	return seq<0,ndim>::all_of([this](int d) -> bool {
		auto& dist(nd[d]);
		if(dist.size() != static_cast<size_t>(np[d] + 1))
			return false;
		if(dist[0] != 0 || dist[np[d]] != n[d])
			return false;
		for(int k = 0; k < np[d]; k++)
			if(dist[k+1] < dist[k])
				return false;
		return true;
	});
}

#if defined(SHARK_PTHREAD_SCHED)

namespace {
	// Auxiliary functions for tdistribution
	template<int d,int ndim>
	INLINE typename enable_if<d < ndim>::type tdistd(const array<int,ndim>& nt, const coords_range<ndim>& l, coords_range<ndim>& t, vector<coords_range<ndim>>& tdist);
	template<int d,int ndim>
	INLINE typename enable_if<d == ndim>::type tdistd(const array<int,ndim>&, const coords_range<ndim>&, coords_range<ndim>& t, vector<coords_range<ndim>>& tdist);

	template<int d,int ndim>
	inline typename enable_if<d < ndim>::type tdistd(const array<int,ndim>& nt, const coords_range<ndim>& l, coords_range<ndim>& t, vector<coords_range<ndim>>& tdist) {
		coord count = l.upper[d] - l.lower[d];
		coord div = count / nt[d];
		coord rem = count % nt[d];

		for(int k = 0; k < nt[d]; k++) {
			t.lower[d] = l.lower[d] + k * div + (k <= rem ? k : rem);
			t.upper[d] = l.lower[d] + (k+1) * div + ((k+1) <= rem ? k+1 : rem);
			tdistd<d+1,ndim>(nt, l, t, tdist);
		}
	}

	template<int d,int ndim>
	inline typename enable_if<d == ndim>::type tdistd(const array<int,ndim>&, const coords_range<ndim>&, coords_range<ndim>& t, vector<coords_range<ndim>>& tdist) {
		tdist.push_back(t);
	}
}

template<int ndim>
vector<coords_range<ndim>> Domain<ndim>::tdistribution() const {
	array<int,ndim> nt = {{}};
	assign_dims(nthrds, ndim, nt.data());
	
	vector<coords_range<ndim>> tdist;
	coords_range<ndim> t;
	tdistd<0,ndim>(nt, local(), t, tdist);

	assert(tdist.size() == static_cast<size_t>(nthrds));

	return tdist;
}
#endif

template<int ndim>
Domain<ndim>::Domain(const Group& group, coords<ndim> n, periods pd, pcoords np)
	: group(group)
	, n(n)
	, pd(pd)
	// Using comma operator to invoke adjustProcs
	, np((adjustProcs(group.nprocs,np), np))
	, nd(distribution(n, np))
	, b(base(np))
#if defined(SHARK_PTHREAD_SCHED)
	, tdist(tdistribution())
#elif defined(SHARK_TBB_SCHED)
	, ap(new tbb::affinity_partitioner())
#endif
{
	assert(b[0] == group.nprocs);
}

template<int ndim>
Domain<ndim>::Domain(const Group& group, coords<ndim> n, periods pd, pcoords np, dists nd)
	: group(group)
	, n(n)
	, pd(pd)
	, np(np)
	, nd(nd)
	, b(base(np))
#if defined(SHARK_PTHREAD_SCHED)
	, tdist(tdistribution())
#elif defined(SHARK_TBB_SCHED)
	, ap(new tbb::affinity_partitioner())
#endif
{
	assert(b[0] == group.nprocs);
	assert(consistentDistribution());
}

template<int ndim>
Domain<ndim>::~Domain() { }

template<int ndim>
void Domain<ndim>::outputDistribution(ostream& out) const {
	out << "Domain " << n << endl;

	for(int k = 0; k < group.nprocs; k++)
		out << "Process " << k << " owns " << local(k) << endl;
}

template<int ndim>
void Domain<ndim>::findi(const vector<coord>& dist, coord i, int& id, coord& off) {
	assert(i >= 0 && i <= dist.back());
	int lo(0), up(static_cast<int>(dist.size()));
	while (lo + 1 < up) {
		int mi((lo + up) / 2);
		if(i >= dist[mi])
			lo = mi;
		else
			up = mi;
	}
	id = lo;
	off = i - dist[lo];
}

template<int ndim>
void Domain<ndim>::find(coords<ndim> i, pcoords& ip, coords<ndim>& off) const {
	for(int d = 0; d < ndim; d++) {
		Domain<ndim>::findi(nd[d], i[d], ip[d], off[d]);
	}
	// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51494
	//seq<0,ndim>::for_each([this,&i,&ip,&off](int d) {
	//	const std::vector<coord>& dist(nd[d]);
	//	findi(dist, i[d], ip[d], off[d]);
	//});
}

template<int ndim>
bool Domain<ndim>::operator==(const Domain<ndim>& other) const {
	if(this == &other)
		return true;

	return equiv(other) && np == other.np && nd == other.nd;
}

template<int ndim>
bool Domain<ndim>::equiv(const Domain<ndim>& other) const {
	return n == other.n && pd == other.pd && group == other.group;
}

template<int ndim>
void Domain<ndim>::sync() const {
	group.sync();
}

template<int ndim>
int Domain<ndim>::shiftd(int d, int disp, int id) const {
	pcoords ip = indexp(id);
	if(pd[d]) {
		ip[d] = (ip[d] + disp) % np[d];
		if(ip[d] < 0)
			ip[d] += np[d];
		return pindex(ip);
	} else {
		ip[d] += disp;
		if(ip[d] < 0 || ip[d] >= np[d])
#if defined(SHARK_MPI_COMM)
			return MPI_PROC_NULL;
#elif defined(SHARK_NO_COMM)
			return -1;
#else
#error "No null process"
#endif
		else
			return pindex(ip);
	}
}

// Set-up instantiations

#define SYMBD(d) template class Domain<d>;
#include "inst_dim"
#undef SYMBD

