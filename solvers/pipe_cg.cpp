
#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <cassert>
#include <shark.hpp>
#include "pipe_cg.hpp"

using namespace std;
using namespace shark;
using namespace shark::types2d;

int main(int argc, char **argv)
{
	Init(&argc, &argv);

	//Default values that might be overridden by options

	int n = 40; //incremented by one internally
	double tol = 1e-6;
	int maxit = 1000;
	int restart = 30;
	int reps = 1;
	int l = 2;
	ostream* ostr = NULL;

	//Process arguments
	int ch;

	while ((ch = getopt(argc, argv, "n:e:m:r:l:t:h:")) != -1)
	{
		switch (ch)
		{
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

		case 'r':
		{
			istringstream iss(optarg);
			iss >> reps;
		}
			break;

		case 'v':
		{
			ostr = world().procid == 0 ? &cout : NULL;
		}
			break;

		case 'h':
		case '?':
		default:
			cerr << "Usage: " << argv[0]
					<< " [-n <size>] [-e <tolerance>] [-m <max-iterations>]  [-r <nr-solves>] [-v <verbose>] [-t <threads>] [-c cg|cg3|cgd|cg3d|p1cg|cgsgmres|mgsgmres|p1gmres]\n"
					<< endl;
			Abort(1);
		}
	}

	if (world().procid == 0)
	{
		cerr << "sched: " << sched << endl;
		cerr << "n: " << n << endl;
		cerr << "maxit: " << maxit << endl;
		cerr << "tol: " << tol << endl;
		cerr << "nprocs: " << world().nprocs << endl;
		cerr << "nthrds: " << nthrds << endl;
		cerr << "repetitions: " << reps << endl;
	}

	SetupThreads();

	coords size = { { n, n } };
	coords gw = { { 1, 1 } };

	Domain d(world(), size);

	if (world().procid == 0)	d.outputDistribution(cerr);

	GlobalArrayD x_exact(d, gw);

	x_exact = constant(d, 1.0);

	GlobalArrayD b(d, gw);

	GlobalArrayD sol(x_exact, false);

	int k;

	double mintime = 1.0e10;

	for (int i = 0; i < reps; i++)
	{
		applyOperator(b, x_exact);

		sol = constant(d, 0.0);

		double starttime = Wtime();

		pipe_cg(sol, b, tol, k, maxit, 0, ostr);

		double totaltime = Wtime() - starttime;

		if (totaltime < mintime)	mintime = totaltime;
	}

	GlobalArrayD error(sol, true);
	error = error - x_exact;
	double errorNorm = norm2(error);

	if (world().procid == 0)
	{
		cerr << "runtime: " << mintime << endl;
		cerr << "err: " << errorNorm << endl;
		cerr << "it: " << k << endl;
	}

	Finalize();

	return 0;
}
