
#include <iostream>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::ndim;

int main(int argc, char* argv[]) {
	Init(&argc, &argv);
	SetupThreads();

	cerr << "Hello world from " << world().procid << " of " << world().nprocs << endl;
	world().sync();
	{
		coords<3> c = {{100,100,100}};
		Domain<3> dom(world(), c);
		if(world().procid == 0) {
			dom.outputDistribution(cerr);
			for(int k = 0; k < world().nprocs; k++) {
				coords<3> c = dom.count(k);
				cerr << "Process " << k << " has data " << dom.hasData(k) << ": " << c << endl;
				typename Domain<3>::pcoords ip = dom.indexp(k);
				cerr << "Process " << k << "(" << dom.pindex(ip) << ")" << " has indices " << ip[0] << "," << ip[1] << "," << ip[2] << endl;

			}
		}

		GlobalArray<3,double> ga(dom);
		{
			Access<3,double> a(ga);
			const GlobalArray<3,double>& gac(ga);
			const Access<3,double> b(gac);
		}
	
	}

	Finalize();
}

