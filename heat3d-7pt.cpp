/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#include <cmath>
#include <iostream>
#include <sstream>
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

void heat(const GlobalArrayD& ga, GlobalArrayD& gb, double nu)
{
	// Sync and make sure ga is complete
	ga.update();

        coords left =  { -1,  0,  0};
        coords right = {  1,  0,  0};
        coords above = {  0, -1,  0};
        coords below = {  0,  1,  0};
        coords  back = {  0,  0, -1};
        coords front = {  0,  0,  1};

        gb = ga + nu * (
                        shift(ga, left) + 
                        shift(ga,right) + 
                        shift(ga,above) +
                        shift(ga,below) +
                        shift(ga, back) +
                        shift(ga,front)
                        - 7 * ga);

}

int heat_loop(int n, GlobalArrayD& ga, GlobalArrayD& gb, int nr, double dt) {
	double nu = dt * n * n * n;
        int k;

        double start = Wtime();
        double end = start;

	for(k = 0; (end-start) < nr; k++)
	{
		if (k % 2 == 0)
			heat(ga, gb, nu);
		else
			heat(gb, ga, nu);

                end = Wtime();
	}

        return k;
}

int main(int argc, char **argv)
{
	Init(&argc, &argv);

	//Default values that might be overridden by options
	int nr = 2; // run for 2seconds
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
					iss >> nr; 
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

	if(world().procid == 0)
	{
		cerr << "sched: " << sched << endl;
		cerr << "n: " << n  << " (" << n*n*n*sizeof(double)/1e9 << "GB)" << endl;
		cerr << "nr: " << nr << " sec" << endl;
		cerr << "block: " << boolalpha << block << noboolalpha << endl;
		cerr << "nprocs: " << world().nprocs << endl;
		cerr << "nthrds: " << nthrds << endl;
	}

	//Declare auxiliary variables and run	
	double dt = 0.25 / n / n;
	SetupThreads();

	{
		const coords size  = {{n,n,n}};
		const coords ghost = {{1,1,1}};

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
		int iter = heat_loop(n, ga, gb, nr, dt);
		double endtime = Wtime();
		double e1 = norm1(ga) / n / n;
		double f1 = norm2(ga) / n / n;

		if(world().procid == 0)
		{
			cerr << "e1: " << e1 << endl;
			cerr << "f1: " << f1 << endl;
			cerr << "t: " << dt * nr << endl;
                        cerr << "iter: " << iter << endl;
			cerr << "runtime: " << endtime - starttime << endl;
                        cerr << "MLUP/s: " << (double)iter*n*n*n/(endtime - starttime)/1e6 << endl;
		}
	}

	Finalize();

	return 0;
}

