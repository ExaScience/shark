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
#include "laplaceOperator.hpp"
#include "cg.hpp"
#include "pipe_cg.hpp"
#include "pipe_chronopoulos_cg.hpp"
#include "pipe_gmres.hpp"

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
	string method = "cg";
	BASIS basis = MONOMIAL;
	int restart = 30;
	int reps = 1;
	int l = 2;
	ostream* ostr = NULL;

	//Process arguments
	int ch;

	while ((ch = getopt(argc, argv, "n:e:m:r:l:t:hvc:")) != -1)
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

		case 'l':
		{
			istringstream iss(optarg);
			iss >> l;
		}
			break;

		case 'v':
		{
			ostr = world().procid == 0 ? &cout : NULL;
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
		if (method.compare("pipe_gmres") == 0 ) cerr << "cg-method (l): " << method << " (" << l << ") "<< endl;
		else
			cerr << "cg-method: " << method << endl;
		cerr << "repetitions: " << reps << endl;
	}

	SetupThreads();
	{
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

			if (method.compare("cg") == 0)
			{
				//	cg<double, GlobalArrayD>(sol, b, tol, k, maxit, ostr);
			}
			else if (method.compare("cg3") == 0)
			{
				//	cg3<double, GlobalArrayD>(sol, b, tol, k, maxit, ostr);
			}
			else if (method.compare("chronopoulos_cg") == 0)
			{
				//	chronopoulos_cg<double, GlobalArrayD>(sol, b, tol, k, maxit, ostr);

#ifdef SHARK_MPI_ASYNC
			}
			else if (method.compare("pipe_cg") == 0)
			{
				pipe_cg<double, GlobalArrayD>(sol, b, tol, k, maxit, l, 0.0, 0.0, ostr);
			}
			else if (method.compare("pipe_chronopoulos_cg") == 0)
			{
				pipe_chronopoulos_cg<double, GlobalArrayD>(sol, b, tol, k, maxit, ostr);
			}
			else if (method.compare("pipe_gmres") == 0)
			{
				pipe_gmres<double, GlobalArrayD>(sol, b, tol, k, maxit, restart, l, basis, 2*M_PI*M_PI, 8*(n+1)*(n+1), ostr);
#endif
			}
			else
			{
				cerr << "Invalid CG method" << endl;
				Abort(1);
			}

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

	}

	Finalize();

	return 0;
}
