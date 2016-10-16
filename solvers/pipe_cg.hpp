


void Solver::pipe_cg()
{

    vector<GlobalArrayD*> r_w;

    r_w.push_back(new GlobalArrayD(x, false));
    r_w.push_back(new GlobalArrayD(x, false));

    GlobalArrayD& r = *r_w[0];
    GlobalArrayD& w = *r_w[1];

    GlobalArrayD z(r, false), x_old(r, false), r_old(r, false), w_old(r, false);
    GlobalArrayD tmp(b, false); // Tmp vector to make add operations idempotent

    for (int i = 0; i < reps; i++)
    {
        start();

        applyOperator(tmp, x, degree);
        r <<= b - tmp;
        applyOperator(tmp, r, degree);
        w <<= tmp;

        double rho = 1.0, rho_old = 1.0, mu, mu_old = 0.0, gamma, gamma_old = 0.0, nu, res;

        int k;
        for (k = 0; k < maxit; k++)
        {
            Future<valarray<double> > future_dots;
            {SHARK_COUNTER("VecTDot"); future_dots = cidot(r, r_w); }

            {SHARK_COUNTER("MatMult"); applyOperator(z, w, degree); }

            valarray<double> dots;
            {SHARK_COUNTER("VecTDot"); { SHARK_COUNTER("wait"); dots = future_dots.wait(); }}

            mu = dots[0];
            nu = dots[1];

            res = sqrt(mu);

            if (tol > .0 && res <= tol) break;
            if (out != NULL) *out << k << "\t" << res << std::endl;

            {
                SHARK_COUNTER("VecAXPY");
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

                if (k > 0) tmp <<= (1 - rho) * w_old + rho * w - rho * gamma * z ;
                else       tmp <<= rho * w - rho * gamma * z;

                swap(w, w_old);
                swap(w, tmp);

                mu_old = mu;
                rho_old = rho;
                gamma_old = gamma;
            }

        }

        stop(k);
    }

    for(auto a : r_w) delete a;
}


