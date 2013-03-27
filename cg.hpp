/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __CG_HPP
#define __CG_HPP

template<typename Vector, typename Scalar, typename Func>
void cg(const Func& Amult, Vector& x, const Vector& b, Scalar tol, int& k, int maxit) {

	Vector r(b, false), p(x, false), w(p, false);

	r = b - Amult(x);
	p = r;

	Scalar rho = norm2(r);

	for(k = 0; k < maxit; k++) {
		if(rho <= tol)
			break;

		w = Amult(p);

		Scalar alpha = rho*rho / dot(p,w);
		x = x + alpha * p;
		r = r - alpha * w;

		Scalar rho_old = rho;
		rho = norm2(r);

		Scalar beta = rho*rho / (rho_old*rho_old);
		p = r + beta * p;
	}

}

#endif
