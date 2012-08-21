
#include <iostream>
#include <valarray>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types3d;

int main(int argc, char* argv[]) {
	Init(&argc, &argv);
	SetupThreads();

	cerr << "Hello world from " << world().procid << " of " << world().nprocs << endl;
	world().sync();
	{
		coords c = {{100,100,100}};
		Domain dom(world(), c);
		if(world().procid == 0) {
			dom.outputDistribution(cerr);
			for(int k = 0; k < world().nprocs; k++) {
				coords c = dom.count(k);
				cerr << "Process " << k << " has data " << dom.hasData(k) << ": " << c << endl;
				typename Domain::pcoords ip = dom.indexp(k);
				cerr << "Process " << k << "(" << dom.pindex(ip) << ")" << " has indices " << ip[0] << "," << ip[1] << "," << ip[2] << endl;

			}
		}

		coords gw = {{2,2,2}};
		GlobalArrayD ga(dom, gw);
		{
			AccessD a(ga);
			const GlobalArrayD& gac(ga);
			const AccessD b(gac);
		}
		ga = constant(dom, 0.0);
		coords_range inner = {gw,c-gw};
		ga.region(inner) = constant(dom, 1.23);
		{
			const AccessD a(ga);
			dom.for_each(dom.total(), [&a, &inner](coords i) {
				if(a(i) != (inner.contains(i) ? 1.23 : 0.0))
					cout << i << ": " << a(i) << endl;
			});
		}
		GlobalArrayD gb(dom);
		gb = 2 * ga + (-ga) + ga - ga;
		gb = gb * gb;
		{
			const AccessD b(gb);
			dom.for_each(inner, [&b](coords i) {
					if(b(i) != 1.23 * 1.23)
						cout << i << ": " << b(i) << endl;
			});
		}
		{
			double start = Wtime();
			double aa = sum(0.0, ga * ga);
			double ab = sum(0.0, ga * gb);
			double bb = sum(0.0, gb * gb);
			double end = Wtime();
			if(world().procid == 0) {
				cerr << "aa: " << aa << " ab: " << ab << " bb: " << bb << " in " << (end - start) << endl;
			}
		}
		{
			double start = Wtime();
			const ndim::vec<3, double> zero = {{}};
			const ndim::vec<3, double> sums = sum(zero, as_vec(ga * ga, ga * gb, gb * gb));
			double end = Wtime();
			if(world().procid == 0) {
				cerr << "aa: " << sums[0] << " ab: " << sums[1] << " bb: " << sums[2] << " in " << (end - start) << endl;
			}
		}
		{
			double start = Wtime();
			valarray<double> v(3);
			{
				const AccessD a(ga);
				const AccessD b(gb);
				v = sum(v, nullary(dom, [&a,&b](coords ii) {
					valarray<double> v(3);
					v[0] = a(ii) * a(ii);
					v[1] = a(ii) * b(ii);
					v[2] = b(ii) * b(ii);
					return v;
				}));
			}
			double end = Wtime();
			if(world().procid == 0) {
				cerr << "aa: " << v[0] << " ab: " << v[1] << " bb: " << v[2] << " in " << (end - start) << endl;
			}
		}
	}

	Finalize();
}

