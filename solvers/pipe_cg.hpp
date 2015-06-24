/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <iostream>
#include <cmath>

using namespace shark;
using namespace std;

#ifdef SHARK_MPI_ASYNC

// This makes it easier to only store the last \ell, or the last
// 2*\ell elements. The index should always be positive!

template <typename T> class ModuloArray
{
	int n; T* a;

public:
	ModuloArray() : a(NULL) {}
	//~ModuloArray() { if (a != NULL) delete a; }
	ModuloArray(int N) : n(N), a(new T[N]) {}
	T& operator[](int i) {return a[i % n];}
};

template <typename T> class ModuloArray2
{
	int n1, n2; ModuloArray<T>* a;

public:

	ModuloArray2(int N1, int N2) : n1(N1), n2(N2), a(new ModuloArray<T>[N1])
	{
		for (int j=0; j<N1; j++) a[j] = ModuloArray<T>(N2);
	}

	//~ModuloArray2() { if (a != NULL) delete a; }
	ModuloArray<T> operator[](int i) {return a[i % n1];}
};

template <typename Scalar>
void recurrenceCoefficients(int l, int a, Scalar sigma[], ModuloArray<Scalar>& rho,
		ModuloArray<Scalar> gamma, Scalar& c, ModuloArray<Scalar>& e, ModuloArray2<Scalar>& d, ModuloArray2<Scalar>& f)
{
	Scalar Ba_1a = (a < l) ? 0.0 : (1.0 - rho[a-l]) / (rho[a-l] * gamma[a-l]);

	Scalar Baa = (a < l) ? sigma[a] : rho[a-l] / (rho[a-l] * gamma[a-l]);

	Scalar Ba1a = (a < l) ? 1.0 : -1.0 / (rho[a-l] * gamma[a-l]);

	e[a] = (a == 0) ? 1.0 : -rho[a-1] * gamma[a-1] * c;

	c = e[a] * Ba1a;

	if (a > 0)
	{
		for (int j=max(0,a-2*l); j<a; j++)
			d[j][a] = -rho[a-1] * gamma[a-1] * f[j][a-1];

		if (a > 1) d[a-2][a] += 1 - rho[a-1];
			d[a-1][a] += rho[a-1];
	}

	for (int j=max(0,a-2*l+1); j<=a; j++)
	{
		f[j][a] = 0.0;

		if (j > 0) f[j][a] -= d[j-1][a] / (rho[j-1] * gamma[j-1]);

		if (j < a) f[j][a] += d[ j ][a] / gamma[j];

		if (j < a-1) f[j][a] += d[j+1][a] * (1.0-rho[j+1])/(rho[j+1] * gamma[j+1]);
	}

	for (int j=max(0,a-2*l+1); j<a; j++)
		f[j][a] -= Baa * d[j][a];

	f[a][a] += Baa;

	if (a > 0)
	{
		for (int j=max(0,a-2*l+1); j<a-1; j++)
		f[j][a] -= e[a] / e[a-1] * Ba_1a * d[j][a-1];
		f[a-1][a] += e[a] / e[a-1] * Ba_1a;
	}
}

/**
 * Implementation of pipelined CG with variable pipelining depth.
 * See my working notes, without them the code will be hard to read.
 */

template <typename Scalar, typename LinearVector>
Scalar pipe_cg(LinearVector& x, LinearVector& b, Scalar tol, int& i, int maxit,
		int l, double lmin=0.0, double lmax=0.0,ostream* out = NULL)
{
	LinearVector x_swap(b, false), w(b, false);
	ModuloArray<LinearVector> r(2*l), z(l+1);

	for (int j=0; j<2*l; j++) r[j] = LinearVector(b, false);

	for (int j=0; j<l+1; j++) z[j] = LinearVector(b, false);

	vector<LinearVector> temp;
	ModuloArray<Future<valarray<Scalar> > > future_dots(l);
	vector<LinearVector*> tempVector;

	applyOperator(r[0], x);

	r[0] = b - r[0];

	z[0] = r[0];

	//-- Chebyshev basis information --

	Scalar sigma[l];

	for (int j=0; j<l; j++)
		sigma[j] = 0.5*(lmin+lmax) + 0.5*(lmax-lmin)*cos(M_PI*(2.0*j + 1.0)/(2.0*l));

	Scalar c = 0.0, res, nu;
	ModuloArray<Scalar> e(l), mu(2*l+1), rho(2*l), gamma(2*l);
	ModuloArray2<Scalar> d(2*l, l), f(2*l, 2), G(3*l, l+1);
	d[0][0] = 0.0;
	rho[0] = 1.0;
	G[0][0] = idot(z[0],z[0]); // This can be overlapped with the first matvec

	for (i=0; i<maxit; i++)
	{
		applyOperator(z[i+1], z[i]);

		if (i < l)
			z[i+1] = z[i+1] - sigma[i] * z[i];

		int a = i - l;

		if (a >= 0)
		{
			recurrenceCoefficients(l, a, sigma, rho, gamma, c, e, d, f);

			{ //--- wait for the dot-products and put them in G
				int k=0;
				valarray<Scalar> dots = future_dots[a+1].wait();
				for (int j=max(0,a-3*l); j<max(0,a-l+2); j++) G[j][a+1] = dots[k++];
				for (int j=max(0,a-l+2); j<a+2; j++) G[j][a+1] = dots[k++];
			}

			for (int j=max(0,a-l+2); j<=a; j++)
			{
				G[j][a+1] = e[j] * G[j][a+1];

				for (int j2=max(j-2*l,0); j2<j; j2++)
					G[j][a+1] += d[j2][j] * G[j2][a+1];
			}

			mu[a] = e[a] * e[a] * G[a][a];

			for (int j=max(a-2*l,0); j<a; j++)
				mu[a] += 2 * e[a] * G[j][a] * d[j][a] + d[j][a] * mu[j] * d[j][a];

			res = sqrt(fabs(mu[a]));

			if (mu[a] <= tol*tol) break;

			if (out != NULL) *out << i << "\t" << res << std::endl;

			nu = c * G[a][a+1] + f[a][a] * mu[a];
			gamma[a] = mu[a] / nu;

			if (a > 0)
				rho[a] = 1.0 / (1.0 - gamma[a] / gamma[a-1] * mu[a] / mu[a-1] / rho[a-1]);

			w = c * z[a+1];

			for (int j=max(0,a-2*l+1); j<=a; j++) w = w + f[j][a] * r[j];

			if (a > 0) r[a+1] = (1.0 - rho[a]) * r[a-1] + rho[a] * r[a] -(rho[a] * gamma[a]) * w;
			else r[a+1] = rho[a] * r[a] - (rho[a] * gamma[a]) * w;

			x_swap = (1.0 - rho[a]) * x_swap + rho[a] * x + rho[a] * gamma[a] * r[a];

			z[i+1] = (1.0 - rho[a]) * z[i-1] + rho[a] * z[i] - rho[a] * gamma[a] * z[i+1];

			swap(x, x_swap);
		}

		tempVector.clear();

		for (int j=max(0,a-2*l); j<max(0,a+2); j++) tempVector.push_back(&r[j]);

		for (int j=max(0,a+2); j<i+2; j++) tempVector.push_back(&z[j]);

		future_dots[i+1] = cidot(z[i+1],tempVector);
	}

	if (out != NULL)
	{
		*out << "# PIPE(" << l << ")-CG: iteration " << i << endl;
		*out << "#       residual norm: " << res << endl;
	}
	return res;
}
#endif

