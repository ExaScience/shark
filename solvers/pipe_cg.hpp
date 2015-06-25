
#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <cassert>
#include <shark.hpp>
#include "laplaceOperator.hpp"

using namespace std;
using namespace shark;
using namespace shark::types2d;

#ifdef SHARK_MPI_ASYNC
double pipe_cg(GlobalArrayD& x, GlobalArrayD& b,double tol, int& k, int maxit, double h, ostream* out = NULL)
{
  vector<GlobalArrayD*> r_w;

  r_w.push_back(new GlobalArrayD(x, false));
  r_w.push_back(new GlobalArrayD(x, false));

  GlobalArrayD& r = *r_w[0];
  GlobalArrayD& w = *r_w[1];

  GlobalArrayD z(r, false), x_old(r, false), r_old(r, false), w_old(r, false);
  GlobalArrayD tmp(r, false); // Tmp vector to make add operations idempotent

  applyOperator(tmp, x);

  r = b - tmp;

  applyOperator(w, r);

  double rho = 1.0, rho_old = 1.0, mu, mu_old = 0.0, gamma, gamma_old = 0.0, nu, res;

  for (k = 0; k < maxit; k++)
  {
    Future<valarray<double> > future_dots = cidot(r, r_w);

    applyOperator(z, w);

    valarray<double> dots = future_dots.wait();

    mu = dots[0];
    nu = dots[1];

    res = sqrt(mu);

    if (res <= tol) break;
    if (out != NULL) *out << k << "\t" << res << std::endl;

    gamma = mu / nu;

    if (k > 0) rho = 1.0 / (1.0 - gamma / gamma_old * mu / mu_old / rho_old);

    if (k > 0) tmp = (1 - rho) * x_old + rho * x + rho * gamma * r ;
    else       tmp = rho * x + rho * gamma * r;

   	swap(x, x_old);
    swap(x, tmp);

    if (k > 0) tmp = (1 - rho) * r_old + rho * r - rho * gamma * w ;
    else       tmp = rho * r - rho * gamma * w;

    swap(r, r_old);
    swap(r, tmp);

    if (k > 0) tmp = (1 - rho) * w_old + rho * w - rho * gamma * z ;
    else       tmp = rho * w - rho * gamma * z;

    swap(w, w_old);
    swap(w, tmp);

    mu_old = mu;
    rho_old = rho;
    gamma_old = gamma;

  }

  if (out != NULL)
  {
    *out << "# P1CG: iteration " << k << endl;
    *out << "#       residual norm: " << res << endl;
  }

  return res;
}

#endif

