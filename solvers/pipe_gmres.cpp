
#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <cassert>
#include <shark.hpp>
#include <iostream>
#include <valarray>

#include "pipe_gmres.hpp"

using namespace std;
using namespace shark;
using namespace shark::types2d;

int main(int argc, char **argv)
{
	Init(&argc, &argv);

	//Default values that might be overridden by options

	int n = 40; //incremented by one internally
	double tol = 1e-06;
	int maxit = 1000;
	int restart = 30;
	int reps = 1;
	int l = 2;
	ostream* ostr = NULL;
	BASIS basis = NEWTON;
	bool block = false;
	//Process arguments
	int ch;

	while ((ch = getopt(argc, argv, "n:e:m:r:l:t:h:b:")) != -1)
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
		case 'l':
		{
			istringstream iss(optarg);
			iss >> l;
		}
			break;
		case 'b':
					block = true;
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

	const array<int,2> pcoords = {{ 0, block ? 0 : 1 }};
	Domain d(world(), size, pcoords);

	if (world().procid == 0)	d.outputDistribution(cerr);

	GlobalArrayD x_exact(d, gw);

	x_exact = constant(d, 1.0);

	GlobalArrayD b(d, gw);

	GlobalArrayD sol(x_exact, false);

	int k;

	double mintime = 1.0e10;

	applyOperator(b, x_exact);

	sol = constant(d, 0.0);

	double starttime = Wtime();

	pipe_gmres(sol, b, tol, k, maxit, restart, l, basis, 2*M_PI*M_PI, 8*(n+1)*(n+1), ostr);

	double totaltime = Wtime() - starttime;

	if (totaltime < mintime)	mintime = totaltime;

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
