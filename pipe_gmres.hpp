/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <iostream>
#include <valarray>

using namespace std;
using namespace shark;

enum BASIS {
	MONOMIAL, NEWTON, CHEBYSHEV
};

#ifdef SHARK_MPI_ASYNC
template <typename Scalar, typename LinearVector>

void ritz(LinearVector&, int, Scalar*)
{

}

template <typename Scalar, typename LinearVector>
Scalar pipe_gmres(LinearVector& x, LinearVector& b,	Scalar tol, int& totit, int maxit, int restart, int l, BASIS basis, double lmin=0.0,
		double lmax=0.0,ostream* out = NULL)
{
	if (restart > maxit) restart = maxit;

	Scalar givens_c[restart];
	Scalar givens_s[restart];
	Scalar y[restart];
	Scalar b_[restart+1];
	Scalar hess[restart+1][restart];
	Scalar G[restart+1][restart];
	Scalar R[restart][restart];

	Future<std::valarray<Scalar> > future_dots[l];

	std::valarray<Scalar> dot_prods;

	Scalar sigma[l];

	if (basis == MONOMIAL)
		for (int k=0; k<l; k++)
			sigma[k] = 0.0;

	if (basis == NEWTON)
		ritz(x, l, sigma);

	if (basis == CHEBYSHEV)
		for (int k=0; k<l; k++)
			sigma[k] = 0.5*(lmin+lmax) + 0.5*(lmax-lmin)*cos(M_PI*(2.0*k+1.0)/(2.0*l));

	vector<LinearVector> Z(restart+1);
	vector<LinearVector> Q(restart+1);

	for (int i=0; i<=restart; i++)
	{
		Z[i] = LinearVector(b, false);
		Q[i] = LinearVector(b, false);
	}

	vector<LinearVector*> tempVector;

	Scalar rho;
	bool no_conv = true;
	totit = 0;

	while (no_conv)
	{
		applyOperator(Z[0], x);

		Z[0] = b - Z[0];

		rho = norm2(Z[0]);

		Z[0] = Z[0] / rho;

		Q[0] = Z[0];

		b_[0] = rho;

		for (int i=1; i<=restart; i++) b_[i] = 0.0;

		for (int j=0; j<=restart; j++)
			for (int i=0; i<restart; i++)
			{
				hess[j][i] = 0.0;
				G[j][i] = 0.0;
			}

		int nrit = restart-2*l-1;

		if (out != NULL) *out << totit << "\t" << rho << endl;

		for (int it=0; it<restart+l; it++)
		{
			totit++;

			int i = it-l;

			if (it < restart) applyOperator(Z[it+1], Z[it]);

			if (i < 0)
			{
				if (it < restart) Z[it+1] = Z[it+1] - sigma[it] * Z[it];
			}
			else
			{
				// wait for the dot products started l iterations ago to complete
				dot_prods = future_dots[i % l].wait();

				for (int j=0; j<=i-l+1; j++) G[j][i] = dot_prods[j];

				for (int j=max(i-l+2,0); j<i+2; j++) G[j][i] = dot_prods[j];

				for (int j=max(1,i+2-l); j<=i; j++)
				{
					double t = 0;

					for (int k=0; k<j; k++) t += G[k][j-1] * G[k][i];

					G[j][i] = (G[j][i] - t) / G[j][j-1];
				}

				double t = 0.0;

				for (int k=0; k<=i; k++) t += G[k][i] * G[k][i];

				if (G[i+1][i] < t)
				{
					nrit = i-1;
					break;
				}

				G[i+1][i] = sqrt(G[i+1][i] - t);

				if (i == 0)
				{
					hess[0][0] = G[0][0] + sigma[0];
					hess[1][0] = G[1][0];
				}
				else if (i < l)
				{
					for (int j=0; j<=i; j++)
					{
						double t = 0.0;

						for (int k=0; k<i; k++) t += hess[j][k] * G[k][i-1];

						hess[j][i] = (G[j][i] + sigma[i] * G[j][i-1] - t) / G[i][i-1];
					}

					hess[i+1][i] = G[i+1][i] / G[i][i-1];
				}
				else
				{
					for (int j=0; j<=i; j++)
					{
						double t1 = 0.0, t2 = 0.0;

						for (int k=0; k<i+2-l; k++) t1 += G[j][l-1+k] * hess[k][i-l];

						for (int k=0; k<=i; k++) t2 += hess[j][k] * G[k][i-1];

						hess[j][i] = (t1 - t2) / G[i][i-1];
					}

					hess[i+1][i] = G[i+1][i] * hess[i+1-l][i-l] / G[i][i-1];
				}

				Q[i+1] = Z[i+1];

				for (int j=0; j<=i; j++)
					Q[i+1] = Q[i+1] - G[j][i] * Q[j];

				Q[i+1] = Q[i+1] / G[i+1][i];

				if (it < restart)
				{
					for (int j=0; j<=i; j++) Z[it+1] = Z[it+1] - hess[j][i] * Z[j+l];

					Z[it+1] = Z[it+1] / hess[i+1][i];
				}
			}

			//USING DOTS

			/*if (it < restart)
			{
			 for (int j=0; j<=i+1; j++) G[j][it] = dot(Z[it+1],Q[j]);
			 for (int j=max(i+2,0); j<=it+1; j++) G[j][it] = dot(Z[it+1],Z[j]);
			}
			*/

			if (it < restart)
			{
				tempVector.clear();

				for (int j=0; j<=i+1; j++) tempVector.emplace_back(&Q[j]);

				for (int j=max(i+2,0); j<=it+1; j++) tempVector.emplace_back(&Z[j]);

				future_dots[it % l] = cidot(Z[it+1],tempVector);
			}

			if (i >= 0)
			{
				R[0][i] = hess[0][i];

				for (int j=1; j<=i; j++)
				{
					double gamma = givens_c[j-1] * R[j-1][i] + givens_s[j-1] * hess[j][i];
					R[j][i] = -givens_s[j-1] * R[j-1][i] + givens_c[j-1] * hess[j][i];
					R[j-1][i] = gamma;
				}

				double delta = sqrt(R[i][i] * R[i][i] + hess[i+1][i] * hess[i+1][i]);
				givens_c[i] = R[i][i] / delta;
				givens_s[i] = hess[i+1][i] / delta;
				R[i][i] = givens_c[i] * R[i][i] + givens_s[i] * hess[i+1][i];
				b_[i+1] = -givens_s[i] * b_[i];
				b_[i] = givens_c[i] * b_[i];
				rho = fabs(b_[i+1]);
			}

			if (out != NULL) *out << totit << "\t" << rho << endl;
			if (i >= 0)

			if ((rho < tol) || (totit >= maxit + l))
			{
				no_conv = false;
				nrit = i;
				break;
			}
		}

		for (int i=0; i<l; i++) future_dots[i].wait(); // wait for the remaining futures

		for (int k=nrit; k>=0; k--)
		{
			y[k] = b_[k];

			for (int i=k+1; i<=nrit; i++)
				y[k] -= R[k][i] * y[i];

			y[k] /= R[k][k];
		}

		for (int i=0; i<=nrit; i++)
			x = x + y[i] * Q[i];

	}

	return rho;
}

#endif
