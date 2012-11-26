/*
 * Copyright (c) 2010-2012, Vrije Universiteit Brussel.
 * All rights reserved.
 */

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
		GlobalArrayD gb(dom), gc(dom);
		gb = 2 * ga + (-ga) + ga - ga;
		gb = gb * gb;
		{
			const AccessD b(gb);
			dom.for_each(inner, [&b](coords i) {
					if(b(i) != 1.23 * 1.23)
						cout << i << ": " << b(i) << endl;
			});
		}
		gc = abs(ga - gb);
		ga = coord_val<2>(dom, 1.0);
		{
			double start = Wtime();
			double aa = sum(0.0, ga * ga);
			double ab = sum(0.0, ga * gb);
			double ac = sum(0.0, ga * gc);
			double end = Wtime();
			if(world().procid == 0) {
				cerr << "aa: " << aa << " ab: " << ab << " ac: " << ac << " in " << (end - start) << endl;
			}
		}
		{
			double start = Wtime();
			const ndim::vec<3, double> zero = {{}};
			const ndim::vec<3, double> sums = sum(zero, as_vec(ga * ga, ga * gb, ga * gc));
			double end = Wtime();
			if(world().procid == 0) {
				cerr << "aa: " << sums[0] << " ab: " << sums[1] << " ac: " << sums[2] << " in " << (end - start) << endl;
			}
		}
		{
			double start = Wtime();
			vector<GlobalArrayD> vg;
			vg.push_back(std::move(ga));
			vg.push_back(std::move(gb));
			vg.push_back(std::move(gc));
			GlobalArrayD& gar(vg[0]);
			valarray<double> v(vg.size());
			{
				const AccessD a(static_cast<const GlobalArrayD&>(gar));
				vector<AccessD> va;
				for(size_t k = 0; k < vg.size(); k++)
					va.push_back(AccessD(static_cast<const GlobalArrayD&>(vg[k])));
				v = dom.sum(dom.total(), v, [&a,&va](valarray<double>& acc, coords ii) {
					for(size_t k = 0; k < va.size(); k++) {
						acc[k] += a(ii) * va[k](ii);
					}
				});
			}
			double end = Wtime();
			if(world().procid == 0) {
				cerr << "aa: " << v[0] << " ab: " << v[1] << " bb: " << v[2] << " in " << (end - start) << endl;
			}
		}
		{
			GlobalArrayP gap(dom);
			gap = nullary(dom, [](coords ii) -> partD {
				vecD v = {{1.25,1.25,1.25}};
				partD p = {{{0,0,0}}, ii.to_vec() * v};
				return p;
			});
			double start = Wtime();
			{
				AccessP ap(gap);
				dom.for_each(dom.total(), [&ap](coords ii) {
					ap(ii).x += ap(ii).v * 1.23;
				});
			}
			double end = Wtime();
			if(world().procid == 0)
				cerr << (end - start) << endl;
		}
	}

	Finalize();
}

