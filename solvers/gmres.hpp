
void Solver::gmres()
{

   if (restart > maxit)  restart = maxit;

    double givens_c[restart];
    double givens_s[restart];
    double y[restart];
    double b_[restart+1];
    double hess[restart+1][restart];
    double rho, gamma;

    vector<GlobalArrayD> V(restart+1);
    vector<GlobalArrayD*> tempVector(restart+1);

    for(int i=0; i<reps; ++i) {
        start();

        for (int i=0; i<=restart; i++)
            V[i] = GlobalArrayD(b, false);

        bool no_conv = true;
        int totit = 0;

        while (no_conv)
        {
            applyOperator(V[0], x, degree);

            V[0] = b - V[0];

            rho = norm2(V[0]);

            V[0] =  V[0] / rho;

            b_[0] = rho;

            for (int i = 1; i <= restart; i++) b_[i] = 0.0;

            for (int j = 0; j <= restart; j++)
                for (int i = 0; i < restart; i++)
                    hess[j][i] = 0.0;

            int nrit = restart-1;

            if (out != NULL) *out << totit << "\t" << rho << std::endl;

            for (int it=0; it<restart; it++)
            {
                totit++;

                applyOperator(V[it+1], V[it], degree);

                if (GStype == Solver::CLASSICAL)
                {
                    tempVector.clear();

                    for (int i=0; i<=restart; i++) tempVector.emplace_back(&V[i]);

                    // do all dot products in one reduce
                    valarray<double> dot_products = cdot(V[it+1], tempVector);

                    for (int k=0; k<=it; k++)
                        hess[k][it] = dot_products[k];

                    for (int k=0; k<=it; k++)
                        V[it+1] = V[it+1] - hess[k][it] * V[k];
                }
                else if (GStype == Solver::MODIFIED)
                {
                    for (int k=0; k<=it; k++)
                    {
                        hess[k][it] = dot(V[it+1],V[k]);
                        V[it+1] = V[it+1] - hess[k][it] * V[k];
                    }
                }

                hess[it+1][it] = norm2(V[it+1]);

                V[it+1] = V[it+1] / hess[it+1][it];

                for (int k=1; k<it+1; k++)
                {
                    gamma = givens_c[k-1]*hess[k-1][it] + givens_s[k-1]*hess[k][it];
                    hess[k][it] = -givens_s[k-1]*hess[k-1][it] + givens_c[k-1]*hess[k][it];
                    hess[k-1][it] = gamma;
                }

                //row 1 c given_s row 0 c given_c
                if( hess[it+1][it] == 0.0)
                {
                    // It is already lower diagonal. Just leave it be....
                    givens_c[it] = 1.0;
                    givens_s[it] = 0.0;
                }
                else if (fabs(hess[it+1][it]) > fabs(hess[it][it]))
                {
                    // The off diagonal entry has a larger
                    // magnitude. Use the ratio of the
                    // diagonal entry over the off diagonal.

                    double tmp = hess[it][it]/hess[it+1][it];
                    givens_s[it] = 1.0/sqrt(1.0+tmp*tmp);
                    givens_c[it] = tmp * givens_s[it];
                }
                else
                {
                    // The off diagonal entry has a smaller
                    // magnitude. Use the ratio of the off
                    // diagonal entry to the diagonal entry.

                    double tmp = hess[it+1][it]/hess[it][it];
                    givens_c[it] = 1.0/sqrt(1.0+tmp*tmp);
                    givens_s[it] = tmp * givens_c[it];
                }

                double tmp = givens_c[it] * hess[it][it] + givens_s[it] * hess[it+1][it];
                hess[it+1][it] = -givens_s[it] * hess[it][it] + givens_c[it] * hess[it+1][it];
                hess[it][it] = tmp;

                tmp = givens_c[it] * b_[it]  + givens_s[it] * b_[it+1] ;
                b_[it+1]  = -givens_s[it] * b_[it]  + givens_s[it] * b_[it+1];
                b_[it]  = tmp;

                rho = fabs(b_[it+1]);

                /*      delta = sqrt(hess[it][it]*hess[it][it] + hess[it+1][it]*hess[it+1][it]);
                        givens_c[it] = hess[it][it] / delta;
                        givens_s[it] = hess[it+1][it] / delta;
                        hess[it][it] = givens_c[it]*hess[it][it] + givens_s[it]*hess[it+1][it];

                        b_[it+1] = -givens_s[it]*b_[it];
                        b_[it] = givens_c[it]*b_[it];
                        rho = fabs(b_[it+1]);
                        */
                if (out != NULL) *out << totit << "\t" << rho << std::endl;

                if ((tol > .0 && rho < tol) || (totit >= maxit))
                {
                    no_conv = false;
                    nrit = it;
                    break;
                }
            }

            for (int k=nrit; k>=0; k--)
            {
                y[k] = b_[k];

                for (int i=k+1; i<=nrit; i++)
                    y[k] -= hess[k][i]*y[i];

                y[k] /= hess[k][k];
            }

            for (int i=0; i<=nrit; i++)
                x = x + y[i] * V[i];
        }
        stop(totit);
    }
}

