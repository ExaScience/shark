
#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <cassert>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types2d;

class LaplaceAcc;

class LaplaceExp {
	friend class LaplaceAcc;
	const double h;
	const GlobalArrayD& gx;
public:
	typedef LaplaceAcc accessor;
	static const int number_of_dimensions = 2;
	LaplaceExp(double h, const GlobalArrayD& gx): h(h), gx(gx) {
		assert(gx.ghost_width()[0] >= 1 && gx.ghost_width()[1] >= 1);
	}
	~LaplaceExp() { }
	const Domain& domain() const {
		return gx.domain();
	}
	coords_range region() const {
		return gx.region();
	}
};

class LaplaceAcc {
	const double h;
	const AccessD x;
	const coord nx;
	const coord ny;
public:
	LaplaceAcc(const LaplaceExp& exp): h(exp.h), x((exp.gx.update(),exp.gx)), nx(exp.domain().n[0]), ny(exp.domain().n[1]) { }
	~LaplaceAcc() { }
	double operator()(coords ii) const {
		coords left;   left[0] = ii[0]-1;  left[1] = ii[1];
		coords right; right[0] = ii[0]+1; right[1] = ii[1];
		coords below; below[0] = ii[0];   below[1] = ii[1]-1;
		coords above; above[0] = ii[0];   above[1] = ii[1]+1;
		//return (4.0 * x(ii) - x(left) - x(below) - x(above) - x(right)) / h / h;
		double r = 4.0 * x(ii);
		if(ii[0] > 0) r -= x(left);
		if(ii[1] > 0) r -= x(below);
		if(ii[1] < ny-1) r -= x(above);
		if(ii[0] < nx-1) r -= x(right);
		r /= h*h;
		return r;
	}
};

LaplaceExp laplace2D(double h, const GlobalArrayD& x) {
	return LaplaceExp(h, x);
}

void cg(double h, GlobalArrayD& x, const GlobalArrayD& b, double tol, int& k, int maxit) {

	GlobalArrayD r(b.domain(), b.ghost_width());
	r = b - laplace2D(h,x);

	GlobalArrayD p(x.domain(), x.ghost_width());
	p = r;

	GlobalArrayD w(p.domain(), p.ghost_width());

	double rho = norm2(r);

	for(k = 0; k < maxit; k++) {
		if(rho <= tol)
			break;

		w = laplace2D(h, p);

		double alpha = rho*rho / dot(p,w);
		x = x + alpha * p;
		r = r - alpha * w;

		double rho_old = rho;
		rho = norm2(r);

		double beta = rho*rho / (rho_old*rho_old);
		p = r + beta * p;
	}

}

int main(int argc, char **argv) {
	Init(&argc, &argv);

	//Default values that might be overridden by options
	int n = 600;	//incremented by one internally
	double tol = 1e0;
	int maxit = 1000;
	string method = "cg";

	//Process arguments
	int ch;
	while((ch = getopt(argc, argv, "n:e:m:t:hc:")) != -1) {
		switch(ch) {
			case 'n':
				{
					istringstream iss(optarg);
					iss >> n;
				}
				break;
			case 'e':
				{
					istringstream iss(optarg);
					iss >> tol; 
				}
				break;
			case 'm':
				{
					istringstream iss(optarg);
					iss >> maxit;
				}
				break;
			case 't':
				{
					istringstream iss(optarg);
					iss >> nthrds;
				}
				break;
			case 'c':
				{
					istringstream iss(optarg);
					iss >> method;
				}
				break;
			case 'h':
			case '?':
			default:
				cerr << "Usage: " << argv[0] << " [-n <size>] [-e <tolerance>] [-m <max-iterations>] [-t <threads>] [-c cg|cgd]\n" << endl;
				Abort(1);
		}
	}

	if(world().procid == 0) {
		cerr << "sched: " << sched << endl;
		cerr << "n: " << n << endl;
		cerr << "maxit: " << maxit << endl;
		cerr << "tol: " << tol << endl;
		cerr << "nprocs: " << world().nprocs << endl;
		cerr << "nthrds: " << nthrds << endl;
		cerr << "cg-method: " << method << endl;
	}

	SetupThreads();
	{
		double h = 1.0 / (n + 1.0);

		coords size = {{n,n}};
		coords gw = {{1,1}};
		Domain d(world(), size);
		if(world().procid == 0)
			d.outputDistribution(cerr);

		GlobalArrayD x_exact(d, gw);
		x_exact = constant(d, 1.0);

		GlobalArrayD b(d);
		b = laplace2D(h, x_exact);

		GlobalArrayD sol(d, gw);
		sol = constant(d, 0.0);

		int k;
		double starttime = Wtime();
		if (method.compare("cg") == 0) {
			cg(h, sol, b, tol, k, maxit);
		} else {
			cerr << "Invalid CG method" << endl;
			Abort(1);
		}
		double endtime = Wtime();

		double err = norm2(sol - x_exact);
		if(world().procid == 0) {
			cerr << "runtime: " << endtime - starttime << endl;
			cerr << "err: " << err << endl;
			cerr << "it: " << k << endl;
		}
	}

	Finalize();

	return 0;
}
