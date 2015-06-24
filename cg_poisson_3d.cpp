/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>
#include <shark.hpp>
#include "cg.hpp"

#define USE_BOUNDS

using namespace std;
using namespace shark;
using namespace shark::types3d;

class LaplaceAcc;

class LaplaceExp {
	friend class LaplaceAcc;
	const GlobalArrayD& gx;
public:
	typedef LaplaceExp storage;
	typedef LaplaceAcc accessor;
	static const int number_of_dimensions = 3;
	LaplaceExp(const GlobalArrayD& gx);
	~LaplaceExp();
	const Domain& domain() const {
		return gx.domain();
	}
	coords_range region() const {
		return gx.region();
	}
};

LaplaceExp::LaplaceExp(const GlobalArrayD& gx): gx(gx) {
	assert(gx.ghost_width()[0] >= 1 && gx.ghost_width()[1] >= 1 && gx.ghost_width()[2] >= 1);
}

LaplaceExp::~LaplaceExp() { }

class LaplaceAcc {
	const AccessD x;
	const coords n;
public:
	LaplaceAcc(const LaplaceExp& exp): x((exp.gx.update(),exp.gx)), n(exp.domain().n) { }
	~LaplaceAcc() { }
	double operator()(coords ii) const {
		coords left;   left[0] = ii[0]-1;  left[1] = ii[1];    left[2] = ii[2];
		coords right; right[0] = ii[0]+1; right[1] = ii[1];   right[2] = ii[2];
		coords below; below[0] = ii[0];   below[1] = ii[1]-1; below[2] = ii[2];
		coords above; above[0] = ii[0];   above[1] = ii[1]+1; above[2] = ii[2];
		coords back;   back[0] = ii[0];    back[1] = ii[1];    back[2] = ii[2]-1;
		coords front; front[0] = ii[0];   front[1] = ii[1];   front[2] = ii[2]+1;

#ifdef USE_BOUNDS
		return 6.0 * x(ii) - x(left) - x(below) - x(above) - x(right) - x(back) - x(front);
#else
		double r = 6.0 * x(ii);
		if(ii[0] > 0) r -= x(left);
		if(ii[0] < n[0]-1) r -= x(right);
		if(ii[1] > 0) r -= x(below);
		if(ii[1] < n[1]-1) r -= x(above);
		if(ii[2] > 0) r -= x(back);
		if(ii[2] < n[2]-1) r -= x(front);
		return r;
#endif
	}
};

int main(int argc, char **argv)
{
	Init(&argc, &argv);

	if(argc != 5)
	{
		cerr << "Usage: " << argv[0] << " <nx> <ny> <nz> <tol>" << endl;
		return 1;
	}
	int nx = stoi(argv[1]);
	int ny = stoi(argv[2]);
	int nz = stoi(argv[3]);
	double tol = stod(argv[4]);

	if(world().procid == 0)
	{
		cerr << "sched: " << sched << endl;
		cerr << "nx: " << nx << endl;
		cerr << "ny: " << ny << endl;
		cerr << "nz: " << nz << endl;
		cerr << "tol: " << tol << endl;
		cerr << "nprocs: " << world().nprocs << endl;
		cerr << "nthrds: " << nthrds << endl;
	}

	SetupThreads();
	{
		auto Amult = [](const GlobalArrayD& x) {
			return LaplaceExp(x);
		};

		coords size = {{nx,ny,nz}};
		coords gw = {{1,1,1}};
#ifdef USE_BOUNDS
		typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};
#else
		typename GlobalArrayD::bounds bd = {{ }};
#endif
		Domain d(world(), size);
		if(world().procid == 0)
			d.outputDistribution(cerr);

		GlobalArrayD x(d, gw, false, bd);
		x = constant(d, 1.0);

		GlobalArrayD b(d);
		b = Amult(x);

		x = constant(d, 0.0);

		int k = 0;
		double starttime = Wtime();
		cg(Amult, x, b, tol, k, 10000);
		double endtime = Wtime();

		if(world().procid == 0)
		{
			cerr << "runtime: " << endtime - starttime << endl;
			cerr << "it: " << k << endl;
		}
	}

	Finalize();

	return 0;
}
