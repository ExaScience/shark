/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <cassert>
#include <shark.hpp>
#include "cg.hpp"

#define USE_BOUNDS

using namespace std;
using namespace shark;
using namespace shark::types2d;

class LaplaceAcc;

class LaplaceExp {
	friend class LaplaceAcc;
	const double h;
	const GlobalArrayD& gx;
public:
	typedef LaplaceExp storage;
	typedef LaplaceAcc accessor;
	static const int number_of_dimensions = 2;
	LaplaceExp(double h, const GlobalArrayD& gx);
	~LaplaceExp();
	const Domain& domain() const {
		return gx.domain();
	}
	coords_range region() const {
		return gx.region();
	}
};

LaplaceExp::LaplaceExp(double h, const GlobalArrayD& gx): h(h), gx(gx) {
	assert(gx.ghost_width()[0] >= 1 && gx.ghost_width()[1] >= 1);
}

LaplaceExp::~LaplaceExp() { }

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
#ifdef USE_BOUNDS
		return (4.0 * x(ii) - x(left) - x(below) - x(above) - x(right)) / h / h;
#else
		double r = 4.0 * x(ii);
		if(ii[0] > 0) r -= x(left);
		if(ii[1] > 0) r -= x(below);
		if(ii[1] < ny-1) r -= x(above);
		if(ii[0] < nx-1) r -= x(right);
		r /= h*h;
		return r;
#endif
	}
};

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
		auto Amult = [h](const GlobalArrayD& x) {
			return LaplaceExp(h, x);
		};

		coords size = {{n,n}};
		coords gw = {{1,1}};
#ifdef USE_BOUNDS
		typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};
#else
		typename GlobalArrayD::bounds bd = {{ }};
#endif
		Domain d(world(), size);
		if(world().procid == 0)
			d.outputDistribution(cerr);

		GlobalArrayD x_exact(d, gw, false, bd);
		x_exact = constant(d, 1.0);

		GlobalArrayD b(d);
		b = Amult(x_exact);

		GlobalArrayD sol(d, gw, false, bd);
		sol = constant(d, 0.0);

		int k = 0;
		double starttime = Wtime();
		if (method.compare("cg") == 0) {
			cg(Amult, sol, b, tol, k, maxit);
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
