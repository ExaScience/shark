/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <iostream>
#include <fstream>
#include<thread>
#include <stdio.h>


using namespace shark;
using namespace std;

#ifdef SHARK_MPI_ASYNC

template <typename Scalar, typename LinearVector> Scalar pipe_chronopoulos_cg(LinearVector& x, LinearVector& b, Scalar tol, int& k, int maxit,
		ostream* out = NULL)
{
	LinearVector r(b, false);

	applyOperator(r, x);

	r = b - r;

	Future<Scalar> f_gamma_old = idot(r,r);

	LinearVector p(r, true);

	LinearVector w(x, false);

	applyOperator(w, r);

	LinearVector s(w, true);

	Scalar gamma_old = f_gamma_old;

	Scalar res = sqrt(gamma_old);

	if (res <= tol) return res;

	vector<const LinearVector* > r_w;

	r_w.emplace_back(&r);

	r_w.emplace_back(&w);

	Future<Scalar> f_wdotr = idot(w,r);

	LinearVector z(s, false);

	applyOperator(z, s);

	LinearVector q(z, true);

	Scalar alpha = gamma_old / f_wdotr;

	Scalar beta = 0.0;

	Future<valarray<Scalar> >  f_dots;

	valarray<Scalar> dots;

	Scalar gamma , delta ;

	for (k = 0; k < maxit; k++)
	{
		if (res <= tol) break;

		if (out != NULL) *out << k << "\t" << res << std::endl;

		z = z * beta + q;
		s = s * beta + w;
		p = p * beta + r;
		x = x + alpha * p;
		r = r - alpha * s;
		w = w - alpha * z;

		if (k < maxit-1)
		{
			f_dots = cidot(r , r_w);

			applyOperator(q, w);

			dots = f_dots.wait();

			gamma = dots[0];
			delta = dots[1];

			beta = gamma / gamma_old;
			alpha = gamma / (delta - beta/alpha * gamma);
			gamma_old = gamma;
			res = sqrt(gamma_old);
		}
	}

	if(out != NULL)
	{
		*out << "# Pipe Chronopoulos CG: iteration " << k << endl;
		*out << "#                  residual norm: " << res << endl;
	}

	return res;
}
#endif

