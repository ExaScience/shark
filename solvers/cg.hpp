/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * All rights reserved.
 */

#ifndef __CG_HPP
#define __CG_HPP

void Solver::onlydot()
{
    double rho;

    for (int i=0; i<reps; ++i)
    {
        start();
        for(int k = 0; k < maxit; k++)
        {
            b = b - x;
            auto f = isum(0.0, b);
            rho = f.wait();
        }
        stop(i);
    }
    std::cout << "rho: " << rho <<std::endl;
}

void Solver::cg() {
    GlobalArrayD r(b, false), p(b, false), w(p, false);
    GlobalArrayD tmp(r, false); // TMP vector to make add operations idempotent

    for (int i=0; i<reps; ++i)
    {
        start();

        applyOperator(tmp, x, degree);
        r = b - tmp;
        p <<= r;

        double rho = norm2(r);
        double alpha;

        int k;
        for(k = 0; k < maxit; k++)
        {
            if(tol > 0.0 && rho <= tol)
                break;

            {SHARK_COUNTER("MatMult"); applyOperator(w, p, degree); }
            {SHARK_COUNTER("VecTDot"); alpha = rho*rho / dot(p,w); }
            {SHARK_COUNTER("VecAXPY"); x = x + alpha * p;}
            {SHARK_COUNTER("VecAXPY"); r = r - alpha * w;}

            double rho_old = rho;
            {SHARK_COUNTER("norm2"); rho = norm2(r);}

            double beta = rho*rho / (rho_old*rho_old);
            {SHARK_COUNTER("VecAYPX"); p <<= r + beta * p; }
        }

        stop(k);
    }
}

#endif
