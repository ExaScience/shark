/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <functional>
#include <getopt.h>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types3d;

void init_pyramid(GlobalArrayD& ga)
{
	const coords_range outer = { ga.domain().total().lower - ga.ghost_width(), ga.domain().total().upper + ga.ghost_width() };
	const vecD one = {{1.0, 1.0, 1.0}};
	const vecD mid = {{0.5, 0.5, 0.5}};

	ga = 1.0 - max_element(abs(coord_vec(ga.domain(), outer, one) - mid)) / 0.5;
}


template<typename S>
struct Heat_3d_7pt {
    double nu;
    double operator()(const typename S::accessor& u, coords ii) const {
        coords  left = ii + coords({{ -1,  0,  0}});
        coords right = ii + coords({{  1,  0,  0}});
        coords above = ii + coords({{  0, -1,  0}});
        coords below = ii + coords({{  0,  1,  0}});
        coords  back = ii + coords({{  0,  0, -1}});
        coords front = ii + coords({{  0,  0,  1}});

        return u(ii) + nu * (u(left) + u(right) + u(above) + u(below) + u(back) + u(front) - 7 * u(ii));
    }
};

template<typename S>
struct Heat_3d_19pt {
    double nu;
    double operator()(const typename S::accessor& u, coords ii) const {
        coords  left = {{ -1,  0,  0}};
        coords right = {{  1,  0,  0}};
        coords above = {{  0, -1,  0}};
        coords below = {{  0,  1,  0}};
        coords  back = {{  0,  0, -1}};
        coords front = {{  0,  0,  1}};

        return u(ii) + nu * (
                u(ii + left*1) + u(ii + right*1) + u(ii + above*1) + u(ii + below*1) + u(ii + back*1) + u(ii + front*1)
              + u(ii + left*2) + u(ii + right*2) + u(ii + above*2) + u(ii + below*2) + u(ii + back*2) + u(ii + front*2)
              + u(ii + left*3) + u(ii + right*3) + u(ii + above*3) + u(ii + below*3) + u(ii + back*3) + u(ii + front*3)
              -  19 * u(ii));
    }
};

void heat_overlap(GlobalArrayD& ga, GlobalArrayD& gb, double nu)
{
    auto op = unary(ga, Heat_3d_19pt<GlobalArrayD>({nu})); 
    auto f = ga.iupdate();
    gb.region(gb.inner()) = op;
    f.wait();
    gb.region(gb.outer()) = op;
}

void heat(GlobalArrayD& ga, GlobalArrayD& gb, double nu)
{
    ga.update();
    gb = unary(ga, Heat_3d_7pt<GlobalArrayD>({nu}));
}

void heat_loop(int n, GlobalArrayD& ga, GlobalArrayD& gb, int iter, double dt) {
	double nu = dt * n * n * n;

	for(int k = 0; k < iter; k+=2)
	{
		heat_overlap(ga, gb, nu);
		heat_overlap(gb, ga, nu);
	}
}

int main(int argc, char **argv)
{
	Init(&argc, &argv);

	//Default values that might be overridden by options
	int iter = -1; 
	int n = 511;	//incremented by one internally
	bool block = false;

	//Process arguments
	int ch;

	while((ch = getopt(argc, argv, "n:i:t:bh")) != -1)
	{
		switch(ch)
		{
			case 'n':
				{
					istringstream iss(optarg);
					iss >> n;
				}
				break;
			case 'i':
				{
					istringstream iss(optarg);
					iss >> iter; 
				}
				break;
			case 't':
				{
					istringstream iss(optarg);
					iss >> nthrds;
				}
				break;
			case 'b':
				block = true;
				break;
			case 'h':
			case '?':
			default:
				cerr << "Usage: " << argv[0] << " [-n <grid-cells>] [-i <iterations>] [-t <threads>]\n" << endl;
				Abort(1);
		}
	}

	// if the number of iterations is not provided we calculate a
	// reasonable number for benchmarking at given size
	if (iter < 0)  iter = 1e7 * world().nprocs / n / n;

	if(world().procid == 0)
	{
		cerr << "sched: " << sched << endl;
		cerr << "n: " << n  << " ( " << n*n*n*sizeof(double)/1e9 << "GB)" << endl;
		cerr << "niter: " << iter << endl;
		cerr << "block: " << boolalpha << block << noboolalpha << endl;
		cerr << "nprocs: " << world().nprocs << endl;
		cerr << "nthrds: " << nthrds << endl;
	}

	//Declare auxiliary variables and run	
	double dt = 0.25 / n / n;
	SetupThreads();

	{
		const coords size  = {{n,n,n}};
		const coords ghost = {{3,3,3}};

		const array<int,3> pcoords = {{ 0, 0, block ? 0 : 1 }};
		Domain d(world(), size, pcoords);

		if(world().procid == 0)
			d.outputDistribution(cerr);

		typename GlobalArrayD::bounds bd = {{ BoundaryD::constant(0.0), BoundaryD::constant(0.0) }};

		GlobalArrayD ga(d, ghost, false, bd);
		GlobalArrayD gb(d, ghost, false, bd);
		world().sync();

		init_pyramid(ga);

		double e0 = norm1(ga) / n / n;
		double f0 = norm2(ga) / n / n;

		if(world().procid == 0)
		{
			cerr << "e0: " << e0 << endl;
			cerr << "f0: " << f0 << endl;
		}

		double starttime = Wtime();
		heat_loop(n, ga, gb, iter, dt);
		double endtime = Wtime();
		double e1 = norm1(ga) / n / n;
		double f1 = norm2(ga) / n / n;

		if(world().procid == 0)
		{
			cerr << "e1: " << e1 << endl;
			cerr << "f1: " << f1 << endl;
			cerr << "t: " << dt * iter << endl;
                        cerr << "iter: " << iter << endl;
			cerr << "runtime: " << endtime - starttime << endl;
                        cerr << "MLUP/s: " << (double)iter*n*n*n/(endtime - starttime)/1e6 << endl;
		}
	}

	Finalize();

	return 0;
}

